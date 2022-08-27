/*
 *  Copyright (c) 2008, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* Author: Brian Gerkey */

#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <memory>

#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>

// Signal handling
#include <signal.h>

#include "amcl/map/map.h"
#include "amcl/pf/pf.h"
#include "amcl/sensors/amcl_odom.h"
#include "amcl/sensors/amcl_laser.h"
#include "portable_utils.hpp"

#include "ros/assert.h"

// roscpp
#include "ros/ros.h"

// Messages that I need
#include "sensor_msgs/LaserScan.h"
#include "geometry_msgs/PoseWithCovarianceStamped.h"
#include "geometry_msgs/PoseArray.h"
#include "geometry_msgs/Pose.h"
#include "geometry_msgs/PoseStamped.h"
#include "nav_msgs/GetMap.h"
#include "nav_msgs/SetMap.h"
#include "std_srvs/Empty.h"

// For transform support
#include "tf2/LinearMath/Transform.h"
#include "tf2/convert.h"
#include "tf2/utils.h"
#include "tf2_geometry_msgs/tf2_geometry_msgs.h"
#include "tf2_ros/buffer.h"
#include "tf2_ros/message_filter.h"
#include "tf2_ros/transform_broadcaster.h"
#include "tf2_ros/transform_listener.h"
#include "message_filters/subscriber.h"

// Dynamic_reconfigure
#include "dynamic_reconfigure/server.h"
#include "amcl/AMCLConfig.h"

// Allows AMCL to run from bag file
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <boost/foreach.hpp>

// For monitoring the estimator
#include <diagnostic_updater/diagnostic_updater.h>

#define NEW_UNIFORM_SAMPLING 1

using namespace amcl;

// Pose hypothesis
typedef struct
{
  // Total weight (weights sum to 1)
  double weight;

  // Mean of pose esimate
  pf_vector_t pf_pose_mean;

  // Covariance of pose estimate
  pf_matrix_t pf_pose_cov;

} amcl_hyp_t;

static double normalize(double z)
{
  return atan2(sin(z),cos(z));
}

static double angle_diff(double a, double b)
{
  double d1, d2;
  a = normalize(a);
  b = normalize(b);
  d1 = a-b;
  d2 = 2*M_PI - fabs(d1);
  if(d1 > 0)
    d2 *= -1.0;
  if(fabs(d1) < fabs(d2))
    return(d1);
  else
    return(d2);
}

static const std::string scan_topic_ = "scan";

/* This function is only useful to have the whole code work
 * with old rosbags that have trailing slashes for their frames
 */
inline std::string stripSlash(const std::string& in)
{
  std::string out = in;
  if ( ( !in.empty() ) && (in[0] == '/') )
    out.erase(0,1);
  return out;
}

class AmclNode
{
  public:
    AmclNode();
    ~AmclNode();

    /**
     * @brief Uses TF and LaserScan messages from bag file to drive AMCL instead
     * @param in_bag_fn input bagfile
     * @param trigger_global_localization whether to trigger global localization
     * before starting to process the bagfile
     */
    void runFromBag(const std::string &in_bag_fn, bool trigger_global_localization = false);

    int process();
    void savePoseToServer();

  private:
    std::shared_ptr<tf2_ros::TransformBroadcaster> tfb_;
    std::shared_ptr<tf2_ros::TransformListener> tfl_;
    std::shared_ptr<tf2_ros::Buffer> tf_;

    bool sent_first_transform_;

    tf2::Transform latest_tf_;
    bool latest_tf_valid_;

    // Pose-generating function used to uniformly distribute particles over
    // the map
    static pf_vector_t uniformPoseGenerator(void* arg);
#if NEW_UNIFORM_SAMPLING
    static std::vector<std::pair<int,int> > free_space_indices;
#endif
    // Callbacks
    bool globalLocalizationCallback(std_srvs::Empty::Request& req,
                                    std_srvs::Empty::Response& res);
    bool nomotionUpdateCallback(std_srvs::Empty::Request& req,
                                    std_srvs::Empty::Response& res);
    bool setMapCallback(nav_msgs::SetMap::Request& req,
                        nav_msgs::SetMap::Response& res);

    void laserReceived(const sensor_msgs::LaserScanConstPtr& laser_scan);
    void initialPoseReceived(const geometry_msgs::PoseWithCovarianceStampedConstPtr& msg);
    void handleInitialPoseMessage(const geometry_msgs::PoseWithCovarianceStamped& msg);
    void mapReceived(const nav_msgs::OccupancyGridConstPtr& msg);

    void handleMapMessage(const nav_msgs::OccupancyGrid& msg);
    void freeMapDependentMemory();
    map_t* convertMap( const nav_msgs::OccupancyGrid& map_msg );
    void updatePoseFromServer();
    void applyInitialPose();

    //parameter for which odom to use
    std::string odom_frame_id_;

    //paramater to store latest odom pose
    geometry_msgs::PoseStamped latest_odom_pose_;

    //parameter for which base to use
    std::string base_frame_id_;
    std::string global_frame_id_;

    bool use_map_topic_;
    bool first_map_only_;

    ros::Duration gui_publish_period;
    ros::Time save_pose_last_time;
    ros::Duration save_pose_period;

    geometry_msgs::PoseWithCovarianceStamped last_published_pose;

    map_t* map_;
    char* mapdata;
    int sx, sy;
    double resolution;

    message_filters::Subscriber<sensor_msgs::LaserScan>* laser_scan_sub_;
    tf2_ros::MessageFilter<sensor_msgs::LaserScan>* laser_scan_filter_;
    ros::Subscriber initial_pose_sub_;
    std::vector< AMCLLaser* > lasers_;
    std::vector< bool > lasers_update_;
    std::map< std::string, int > frame_to_laser_;

    // 粒子滤波相关变量
    pf_t *pf_;
    double pf_err_, pf_z_;
    bool pf_init_;
    pf_vector_t pf_odom_pose_;
    double d_thresh_, a_thresh_;
    int resample_interval_;
    int resample_count_;
    double laser_min_range_;
    double laser_max_range_;

    // 是否强制更新
    bool m_force_update;  // used to temporarily let amcl update samples even when no motion occurs...

    AMCLOdom* odom_;
    AMCLLaser* laser_;

    ros::Duration cloud_pub_interval;
    ros::Time last_cloud_pub_time;

    // For slowing play-back when reading directly from a bag file
    ros::WallDuration bag_scan_period_;

    void requestMap();

    // Helper to get odometric pose from transform system
    bool getOdomPose(geometry_msgs::PoseStamped& pose,
                     double& x, double& y, double& yaw,
                     const ros::Time& t, const std::string& f);

    //time for tolerance on the published transform,
    //basically defines how long a map->odom transform is good for
    ros::Duration transform_tolerance_;

    ros::NodeHandle nh_;
    ros::NodeHandle private_nh_;
    ros::Publisher pose_pub_;
    ros::Publisher particlecloud_pub_;
    ros::ServiceServer global_loc_srv_;
    ros::ServiceServer nomotion_update_srv_; //to let amcl update samples without requiring motion
    ros::ServiceServer set_map_srv_;
    ros::Subscriber initial_pose_sub_old_;
    ros::Subscriber map_sub_;

    diagnostic_updater::Updater diagnosic_updater_;
    void standardDeviationDiagnostics(diagnostic_updater::DiagnosticStatusWrapper& diagnostic_status);
    double std_warn_level_x_;
    double std_warn_level_y_;
    double std_warn_level_yaw_;

    amcl_hyp_t* initial_pose_hyp_;
    bool first_map_received_;
    bool first_reconfigure_call_;

    boost::recursive_mutex configuration_mutex_;
    dynamic_reconfigure::Server<amcl::AMCLConfig> *dsrv_;
    amcl::AMCLConfig default_config_;
    ros::Timer check_laser_timer_;

    int max_beams_, min_particles_, max_particles_;
    double alpha1_, alpha2_, alpha3_, alpha4_, alpha5_;
    double alpha_slow_, alpha_fast_;
    double z_hit_, z_short_, z_max_, z_rand_, sigma_hit_, lambda_short_;
  //beam skip related params
    bool do_beamskip_;
    double beam_skip_distance_, beam_skip_threshold_, beam_skip_error_threshold_;
    double laser_likelihood_max_dist_;
    odom_model_t odom_model_type_;
    double init_pose_[3];
    double init_cov_[3];
    laser_model_t laser_model_type_;
    bool tf_broadcast_;
    bool selective_resampling_;

    void reconfigureCB(amcl::AMCLConfig &config, uint32_t level);

    ros::Time last_laser_received_ts_;
    ros::Duration laser_check_interval_;
    void checkLaserReceived(const ros::TimerEvent& event);
};

#if NEW_UNIFORM_SAMPLING
std::vector<std::pair<int,int> > AmclNode::free_space_indices;
#endif

#define USAGE "USAGE: amcl"

boost::shared_ptr<AmclNode> amcl_node_ptr;

void sigintHandler(int sig)
{
  // Save latest pose as we're shutting down.
  amcl_node_ptr->savePoseToServer();
  ros::shutdown();
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "amcl");
  ros::NodeHandle nh;

  //SIGINT表示默认中断处理，即按Crtl+C
  signal(SIGINT, sigintHandler);

  // Make our node available to sigintHandler
  amcl_node_ptr.reset(new AmclNode());

  if (argc == 1)
  {
    // run using ROS input
    ros::spin();
  }
  else if ((argc >= 3) && (std::string(argv[1]) == "--run-from-bag"))
  {
    if (argc == 3)
    {
      amcl_node_ptr->runFromBag(argv[2]);
    }
    else if ((argc == 4) && (std::string(argv[3]) == "--global-localization"))
    {
      amcl_node_ptr->runFromBag(argv[2], true);
    }
  }

  // Without this, our boost locks are not shut down nicely
  amcl_node_ptr.reset();

  // To quote Morgan, Hooray!
  return(0);
}

AmclNode::AmclNode() :
        sent_first_transform_(false),    // bool, 是否发送过坐标变换
        latest_tf_valid_(false),                // bool, 坐标变换是否可用
        map_(NULL),                                  // map_t*, 地图对象指针
        pf_(NULL),                                        // pf_t*, 粒子滤波器对象指针
        resample_count_(0),                   // int, 重采样计数
        odom_(NULL),                               // AMCLOdom*, 里程计对象
        laser_(NULL),                                 // AMCLLaser*, 激光雷达对象
	      private_nh_("~"),                        // ros::NodeHandle, 局部ROS句柄
        initial_pose_hyp_(NULL),         // amcl_hyp_t*, 初始位姿假设
        first_map_received_(false),      // bool, 是否接收到地图
        first_reconfigure_call_(true)   // 初次调用参数重配置 
{
  boost::recursive_mutex::scoped_lock l(configuration_mutex_);

  // amcl节点相关参数
  #pragma region ParmaInit
  private_nh_.param("use_map_topic", use_map_topic_, false);
  private_nh_.param("first_map_only", first_map_only_, false);
  double tmp;
  private_nh_.param("gui_publish_rate", tmp, -1.0);
  gui_publish_period = ros::Duration(1.0/tmp);
  private_nh_.param("save_pose_rate", tmp, 0.5);
  save_pose_period = ros::Duration(1.0/tmp);
  private_nh_.param("update_min_d", d_thresh_, 0.2);
  private_nh_.param("update_min_a", a_thresh_, M_PI/6.0);
  private_nh_.param("odom_frame_id", odom_frame_id_, std::string("odom"));
  private_nh_.param("base_frame_id", base_frame_id_, std::string("base_link"));
  private_nh_.param("global_frame_id", global_frame_id_, std::string("map"));
  private_nh_.param("resample_interval", resample_interval_, 2);
  private_nh_.param("tf_broadcast", tf_broadcast_, true);
  double tmp_tol;
  private_nh_.param("transform_tolerance", tmp_tol, 0.1);

  //运动模型相关参数
  private_nh_.param("odom_alpha1", alpha1_, 0.2);
  private_nh_.param("odom_alpha2", alpha2_, 0.2);
  private_nh_.param("odom_alpha3", alpha3_, 0.2);
  private_nh_.param("odom_alpha4", alpha4_, 0.2);
  private_nh_.param("odom_alpha5", alpha5_, 0.2);

  //雷达测量模型相关参数
  //1.雷达数据本身的参数配置
  private_nh_.param("laser_min_range", laser_min_range_, -1.0);     //雷达最小距离
  private_nh_.param("laser_max_range", laser_max_range_, -1.0);   //雷达最大距离
  private_nh_.param("laser_max_beams", max_beams_, 30);             //每次扫描要使用多少均匀间隔的光束
  //2.雷达beam模型的相关参数
  private_nh_.param("laser_z_hit", z_hit_, 0.95);
  private_nh_.param("laser_z_short", z_short_, 0.1);
  private_nh_.param("laser_z_max", z_max_, 0.05);
  private_nh_.param("laser_z_rand", z_rand_, 0.05);
  private_nh_.param("laser_sigma_hit", sigma_hit_, 0.2);
  private_nh_.param("laser_lambda_short", lambda_short_, 0.1);
  //3.雷达likelihood_field模型的相关参数
  private_nh_.param("laser_likelihood_max_dist", laser_likelihood_max_dist_, 2.0);
  //4.雷达likelihood_field_prob模型的相关参数
  private_nh_.param("do_beamskip", do_beamskip_, false);
  private_nh_.param("beam_skip_distance", beam_skip_distance_, 0.5);
  private_nh_.param("beam_skip_threshold", beam_skip_threshold_, 0.3);
  if (private_nh_.hasParam("beam_skip_error_threshold_"))
  {
    private_nh_.param("beam_skip_error_threshold_", beam_skip_error_threshold_);
  }
  else
  {
    private_nh_.param("beam_skip_error_threshold", beam_skip_error_threshold_, 0.9);
  }

  //粒子滤波器相关参数
  //1.粒子滤波器本身的相关参数
  private_nh_.param("min_particles", min_particles_, 100);
  private_nh_.param("max_particles", max_particles_, 5000);
  private_nh_.param("selective_resampling", selective_resampling_, false);
  //2.Augmented_MCL相关参数
  private_nh_.param("recovery_alpha_slow", alpha_slow_, 0.001);
  private_nh_.param("recovery_alpha_fast", alpha_fast_, 0.1);
  //3.KLD_MCL相关参数
  private_nh_.param("kld_z", pf_z_, 0.99);
  private_nh_.param("kld_err", pf_err_, 0.01);

  // 用于诊断
  private_nh_.param("std_warn_level_x", std_warn_level_x_, 0.2);
  private_nh_.param("std_warn_level_y", std_warn_level_y_, 0.2);
  private_nh_.param("std_warn_level_yaw", std_warn_level_yaw_, 0.1);

  transform_tolerance_.fromSec(tmp_tol);
  #pragma endregion 

  // step 2: 选择测量模型，参考《概率机器人》第6章，默认likelihood_field参考6.4节的似然域模型
  #pragma region LaserModel
  std::string tmp_model_type;
  private_nh_.param("laser_model_type", tmp_model_type, std::string("likelihood_field"));
  if(tmp_model_type == "beam")
    laser_model_type_ = LASER_MODEL_BEAM;
  else if(tmp_model_type == "likelihood_field")
    laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD;
  else if(tmp_model_type == "likelihood_field_prob"){
    laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD_PROB;
  }
  else
  {
    ROS_WARN("Unknown laser model type \"%s\"; defaulting to likelihood_field model",
             tmp_model_type.c_str());
    laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD;
  }
  #pragma endregion

  // step 3: 选择机器人运动学模型，默认是差分机器人，参考《概率机器人》5.4节里程计运动模型
  #pragma region OdomModel
  private_nh_.param("odom_model_type", tmp_model_type, std::string("diff"));
  if(tmp_model_type == "diff")
    odom_model_type_ = ODOM_MODEL_DIFF;
  else if(tmp_model_type == "omni")
    odom_model_type_ = ODOM_MODEL_OMNI;
  else if(tmp_model_type == "diff-corrected")
    odom_model_type_ = ODOM_MODEL_DIFF_CORRECTED;
  else if(tmp_model_type == "omni-corrected")
    odom_model_type_ = ODOM_MODEL_OMNI_CORRECTED;
  else
  {
    ROS_WARN("Unknown odom model type \"%s\"; defaulting to diff model",
             tmp_model_type.c_str());
    odom_model_type_ = ODOM_MODEL_DIFF;
  }
  #pragma endregion

  {
    double bag_scan_period;
    private_nh_.param("bag_scan_period", bag_scan_period, -1.0);
    bag_scan_period_.fromSec(bag_scan_period);
  }

  odom_frame_id_ = stripSlash(odom_frame_id_);
  base_frame_id_ = stripSlash(base_frame_id_);
  global_frame_id_ = stripSlash(global_frame_id_);

  // step 3: 从参数服务器中得到初始位姿和初始协方差
  updatePoseFromServer();

  cloud_pub_interval.fromSec(1.0);
  tfb_.reset(new tf2_ros::TransformBroadcaster());   //一个坐标变换的广播器
  tf_.reset(new tf2_ros::Buffer());                                       //一个坐标变换的历史缓存
  tfl_.reset(new tf2_ros::TransformListener(*tf_));    //一个坐标变换的监听器

  // step 4 : 定义话题,订阅和发布
  #pragma region TopicAndSeriver
  pose_pub_ = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("amcl_pose", 2, true);
  particlecloud_pub_ = nh_.advertise<geometry_msgs::PoseArray>("particlecloud", 2, true);
  //服务global_localization用于获取机器人的全局定位
  global_loc_srv_ = nh_.advertiseService("global_localization", &AmclNode::globalLocalizationCallback,this);
  //服务request_nomotion_update用于手动的触发粒子更新并发布新的位姿估计
  nomotion_update_srv_= nh_.advertiseService("request_nomotion_update", &AmclNode::nomotionUpdateCallback, this);
  //服务set_map用于设定机器人位姿和地图信息。
  set_map_srv_= nh_.advertiseService("set_map", &AmclNode::setMapCallback, this);

  laser_scan_sub_ = new message_filters::Subscriber<sensor_msgs::LaserScan>(nh_, scan_topic_, 100);
  //主要是同时监听激光扫描消息和里程计坐标变换，同步两者的输出
  laser_scan_filter_ =new tf2_ros::MessageFilter<sensor_msgs::LaserScan>(*laser_scan_sub_,
                                                             *tf_,
                                                             odom_frame_id_,
                                                             100,
                                                             nh_);
  laser_scan_filter_->registerCallback(boost::bind(&AmclNode::laserReceived,this, _1));
  initial_pose_sub_ = nh_.subscribe("initialpose", 2, &AmclNode::initialPoseReceived, this);

  // 要么订阅要么通过服务获得地图
  if(use_map_topic_) 
  {
    map_sub_ = nh_.subscribe("map", 1, &AmclNode::mapReceived, this);
    ROS_INFO("Subscribed to map topic.");
  }
   else
  {
    requestMap();
  }
  m_force_update = false;
 #pragma endregion 

  // step 5: 参数动态配置
  #pragma region Other
  dsrv_ = new dynamic_reconfigure::Server<amcl::AMCLConfig>(ros::NodeHandle("~"));
  dynamic_reconfigure::Server<amcl::AMCLConfig>::CallbackType cb = boost::bind(&AmclNode::reconfigureCB, this, _1, _2);
  dsrv_->setCallback(cb);

  //定义一个计时器用于每隔15s检查一次激光雷达的接收数据，如果期间没有收到新的数据给出警告
  laser_check_interval_ = ros::Duration(15.0);
  check_laser_timer_ = nh_.createTimer(laser_check_interval_,boost::bind(&AmclNode::checkLaserReceived, this, _1));

  //最后将粒子滤波器的标准差注册到诊断信息中
  diagnosic_updater_.setHardwareID("None");
  diagnosic_updater_.add("Standard deviation", this, &AmclNode::standardDeviationDiagnostics);
  #pragma endregion
}

#pragma region MapAndPose
// 请求服务static_server提供map，然后调用handleMapMessage处理地图信息
void AmclNode::requestMap()
{
  boost::recursive_mutex::scoped_lock ml(configuration_mutex_);

  // RPC采用客户端/服务端的模式，通过request-response消息模式实现地图获取
  nav_msgs::GetMap::Request  req;
  nav_msgs::GetMap::Response resp;
  ROS_INFO("Requesting the map...");
  // 请求static_map服务，直到成功，该服务在map_server包的map_server节点中
  while(!ros::service::call("static_map", req, resp))
  {
    ROS_WARN("Request for map failed; trying again...");
    ros::Duration d(0.5);
    d.sleep();
  }
  handleMapMessage( resp.map );
}

//地图数据topic
void AmclNode::mapReceived(const nav_msgs::OccupancyGridConstPtr& msg)
{
  //first_map_only 参数服务器，初始化为false
  //first_map_received_  AMCLNode参数，初始划为false，在本段代码后，将其修改为true，表示该函数只运行一次
  if( first_map_only_ && first_map_received_ ) {
    return;
  }

  handleMapMessage( *msg );

  first_map_received_ = true;
}

//服务set_map的回调函数
bool AmclNode::setMapCallback(nav_msgs::SetMap::Request& req,nav_msgs::SetMap::Response& res)
{
  handleMapMessage(req.map);
  handleInitialPoseMessage(req.initial_pose);
  res.success = true;
  return true;
}

// 处理地图信息，地图转换
void AmclNode::handleMapMessage(const nav_msgs::OccupancyGrid& msg)
{
  // step 1 : 加锁保护
  boost::recursive_mutex::scoped_lock cfl(configuration_mutex_);

  ROS_INFO("Received a %d X %d map @ %.3f m/pix\n",msg.info.width,msg.info.height,msg.info.resolution);

  // step 2： msg的frame_id 需要与 global_frame_id_ 一致
  if(msg.header.frame_id != global_frame_id_)
    ROS_WARN("Frame_id of map received:'%s' doesn't match global_frame_id:'%s'. This could cause issues with reading published topics",
             msg.header.frame_id.c_str(),
             global_frame_id_.c_str());

  // step 3： 释放 map_, pf_, odom_, laser_的空间
  freeMapDependentMemory();
  // Clear queued laser objects because they hold pointers to the existing
  lasers_.clear();
  lasers_update_.clear();
  frame_to_laser_.clear();

  // step 4: 转换成标准地图 0 转换成-1(自由空间);100转换成 +1（障碍物）; 其他的 转换成0（未知）
  map_ = convertMap(msg);

// step 5: 将自由空间点的坐标保存下来
#if NEW_UNIFORM_SAMPLING
  // 自由空间的索引
  free_space_indices.resize(0);
  for(int i = 0; i < map_->size_x; i++)
    for(int j = 0; j < map_->size_y; j++)
      if(map_->cells[MAP_INDEX(map_,i,j)].occ_state == -1)
        free_space_indices.push_back(std::make_pair(i,j));
#endif
  // step 6: 创建粒子滤波器
  pf_ = pf_alloc(min_particles_, max_particles_,
                 alpha_slow_, alpha_fast_,
                 (pf_init_model_fn_t)AmclNode::uniformPoseGenerator,
                 (void *)map_);
  pf_set_selective_resampling(pf_, selective_resampling_);
  pf_->pop_err = pf_err_;
  pf_->pop_z = pf_z_;

  // step 7: 初始化粒子滤波器， 从参数服务器获取初始位姿及方差放到pf中
  updatePoseFromServer();
  pf_vector_t pf_init_pose_mean = pf_vector_zero();
  pf_init_pose_mean.v[0] = init_pose_[0];
  pf_init_pose_mean.v[1] = init_pose_[1];
  pf_init_pose_mean.v[2] = init_pose_[2];
  pf_matrix_t pf_init_pose_cov = pf_matrix_zero();
  pf_init_pose_cov.m[0][0] = init_cov_[0];
  pf_init_pose_cov.m[1][1] = init_cov_[1];
  pf_init_pose_cov.m[2][2] = init_cov_[2];
  pf_init(pf_, pf_init_pose_mean, pf_init_pose_cov);
  pf_init_ = false;

  // step 8: 实例化传感器（里程计和激光雷达）
  // Odometry 里程计
  delete odom_;
  odom_ = new AMCLOdom();
  ROS_ASSERT(odom_);
  odom_->SetModel( odom_model_type_, alpha1_, alpha2_, alpha3_, alpha4_, alpha5_ );

  // Laser 激光雷达
  delete laser_;
  laser_ = new AMCLLaser(max_beams_, map_);
  ROS_ASSERT(laser_);
  if(laser_model_type_ == LASER_MODEL_BEAM)
    laser_->SetModelBeam(z_hit_, z_short_, z_max_, z_rand_,
                         sigma_hit_, lambda_short_, 0.0);
  else if(laser_model_type_ == LASER_MODEL_LIKELIHOOD_FIELD_PROB){
    ROS_INFO("Initializing likelihood field model; this can take some time on large maps...");
    laser_->SetModelLikelihoodFieldProb(z_hit_, z_rand_, sigma_hit_,
					laser_likelihood_max_dist_,
					do_beamskip_, beam_skip_distance_,
					beam_skip_threshold_, beam_skip_error_threshold_);
    ROS_INFO("Done initializing likelihood field model.");
  }
  else
  {
    ROS_INFO("Initializing likelihood field model; this can take some time on large maps...");
    laser_->SetModelLikelihoodField(z_hit_, z_rand_, sigma_hit_,
                                    laser_likelihood_max_dist_);
    ROS_INFO("Done initializing likelihood field model.");
  }

  // step 9: 确定初始位姿，最后才做这一步为了防止初始位姿比地图信息早到
  //在接收到地图数据之前，有可能已经接收到了机器人的初始位姿，这里调用函数applyInitialPose处理这一情况
  applyInitialPose();

}

//释放空间
void AmclNode::freeMapDependentMemory()
{
  if( map_ != NULL ) {
    map_free( map_ );
    map_ = NULL;
  }
  if( pf_ != NULL ) {
    pf_free( pf_ );
    pf_ = NULL;
  }
  delete odom_;
  odom_ = NULL;
  delete laser_;
  laser_ = NULL;
}

//将ros中地图格式转换成amcl的地图格式
map_t* AmclNode::convertMap( const nav_msgs::OccupancyGrid& map_msg )
{
  map_t* map = map_alloc();
  ROS_ASSERT(map);

  map->size_x = map_msg.info.width;
  map->size_y = map_msg.info.height;
  map->scale = map_msg.info.resolution;
  map->origin_x = map_msg.info.origin.position.x + (map->size_x / 2) * map->scale;
  map->origin_y = map_msg.info.origin.position.y + (map->size_y / 2) * map->scale;
  // Convert to player format
  map->cells = (map_cell_t*)malloc(sizeof(map_cell_t)*map->size_x*map->size_y);
  ROS_ASSERT(map->cells);
  //地图 0 转换成-1(自由空间);100转换成 +1（障碍物）; 其他的转换成0（未知）
  for(int i=0;i<map->size_x * map->size_y;i++)
  {
    if(map_msg.data[i] == 0)
      map->cells[i].occ_state = -1;
    else if(map_msg.data[i] == 100)
      map->cells[i].occ_state = +1;
    else
      map->cells[i].occ_state = 0;
  }

  return map;
}

//初始化位姿的回调函数
void AmclNode::initialPoseReceived(const geometry_msgs::PoseWithCovarianceStampedConstPtr& msg)
{
  handleInitialPoseMessage(*msg);
}

//处理初始化位姿的函数
void AmclNode::handleInitialPoseMessage(const geometry_msgs::PoseWithCovarianceStamped& msg)
{
  boost::recursive_mutex::scoped_lock prl(configuration_mutex_);
  if(msg.header.frame_id == "")
  {
    // This should be removed at some point
    ROS_WARN("Received initial pose with empty frame_id.  You should always supply a frame_id.");
  }
  // We only accept initial pose estimates in the global frame, #5148.
  else if(stripSlash(msg.header.frame_id) != global_frame_id_)
  {
    ROS_WARN("Ignoring initial pose in frame \"%s\"; initial poses must be in the global frame, \"%s\"",
             stripSlash(msg.header.frame_id).c_str(),
             global_frame_id_.c_str());
    return;
  }

  // In case the client sent us a pose estimate in the past, integrate the
  // intervening odometric change.
  geometry_msgs::TransformStamped tx_odom;
  try
  {
    ros::Time now = ros::Time::now();
    /*
    实际上，完整的lookupTransform共包括6个入口参数，说明如下：
    lookupTransform(string& target_frame, ros::Time& target_time
            string& source_frame, ros::Time& source_time
            string& fixed_frame, ros::Time& timeout)

       target_frame:目标坐标系名称
       target_time：目标坐标系的时戳
       source_frame:原坐标系名称
       source_time:原坐标系时戳
       fixed_frame:固定坐标系名称
       timeout:等待时间
       ros::Time::now()：该参数表示从tf_buffer中获取和当前时戳一致的tf，由于传输延迟的存在，期望的tf很可能并不存在；
       ros::Time(0)：该参数表示从tf_buffer中获取时戳最新的tf，期望的tf一定存在，但由于传输延迟的存在，该tf并不一定是当前时刻的tf；7
      在默认调用中，target_time和source_time默认一致，固定坐标系默认为世界坐标系"/world"。
    */
    tx_odom = tf_->lookupTransform(base_frame_id_, msg.header.stamp,
                                   base_frame_id_, ros::Time::now(),
                                   odom_frame_id_, ros::Duration(0.5));
  }
  catch(tf2::TransformException e)
  {
    // If we've never sent a transform, then this is normal, because the
    // global_frame_id_ frame doesn't exist.  We only care about in-time
    // transformation for on-the-move pose-setting, so ignoring this
    // startup condition doesn't really cost us anything.

    //sent_first_transform_ AMCLNode 初始化为false,在laserReceived函数中正确广播map到odom的tf树将其置为true
    if(sent_first_transform_)
      ROS_WARN("Failed to transform initial pose in time (%s)", e.what());
    tf2::convert(tf2::Transform::getIdentity(), tx_odom.transform);
  }

  tf2::Transform tx_odom_tf2;
  tf2::convert(tx_odom.transform, tx_odom_tf2);
  tf2::Transform pose_old, pose_new;
  tf2::convert(msg.pose.pose, pose_old);
  pose_new = pose_old * tx_odom_tf2;

  // Transform into the global frame
  ROS_INFO("Setting pose (%.6f): %.3f %.3f %.3f",
           ros::Time::now().toSec(),
           pose_new.getOrigin().x(),
           pose_new.getOrigin().y(),
           tf2::getYaw(pose_new.getRotation()));
  // Re-initialize the filter
  pf_vector_t pf_init_pose_mean = pf_vector_zero();
  pf_init_pose_mean.v[0] = pose_new.getOrigin().x();
  pf_init_pose_mean.v[1] = pose_new.getOrigin().y();
  pf_init_pose_mean.v[2] = tf2::getYaw(pose_new.getRotation());
  pf_matrix_t pf_init_pose_cov = pf_matrix_zero();
  // Copy in the covariance, converting from 6-D to 3-D
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<2; j++)
    {
      pf_init_pose_cov.m[i][j] = msg.pose.covariance[6*i+j];
    }
  }
  pf_init_pose_cov.m[2][2] = msg.pose.covariance[6*5+5];

  delete initial_pose_hyp_;
  initial_pose_hyp_ = new amcl_hyp_t();
  initial_pose_hyp_->pf_pose_mean = pf_init_pose_mean;
  initial_pose_hyp_->pf_pose_cov = pf_init_pose_cov;
  applyInitialPose();
}

//通过调用粒子滤波进行初始化位姿操作
void AmclNode::applyInitialPose()
{
  boost::recursive_mutex::scoped_lock cfl(configuration_mutex_);
  if( initial_pose_hyp_ != NULL && map_ != NULL )
   {
      pf_init(pf_, initial_pose_hyp_->pf_pose_mean, initial_pose_hyp_->pf_pose_cov);
      pf_init_ = false;

      delete initial_pose_hyp_;
      initial_pose_hyp_ = NULL;
  }
}

//通过参数服务器获取机器人上一时刻的位姿
void AmclNode::updatePoseFromServer()
{
  init_pose_[0] = 0.0;
  init_pose_[1] = 0.0;
  init_pose_[2] = 0.0;
  init_cov_[0] = 0.5 * 0.5;
  init_cov_[1] = 0.5 * 0.5;
  init_cov_[2] = (M_PI/12.0) * (M_PI/12.0);
  // Check for NAN on input from param server, #5239
  double tmp_pos;
  private_nh_.param("initial_pose_x", tmp_pos, init_pose_[0]);
  if(!std::isnan(tmp_pos))
    init_pose_[0] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial pose X position");
  private_nh_.param("initial_pose_y", tmp_pos, init_pose_[1]);
  if(!std::isnan(tmp_pos))
    init_pose_[1] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial pose Y position");
  private_nh_.param("initial_pose_a", tmp_pos, init_pose_[2]);
  if(!std::isnan(tmp_pos))
    init_pose_[2] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial pose Yaw");
  private_nh_.param("initial_cov_xx", tmp_pos, init_cov_[0]);
  if(!std::isnan(tmp_pos))
    init_cov_[0] =tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial covariance XX");
  private_nh_.param("initial_cov_yy", tmp_pos, init_cov_[1]);
  if(!std::isnan(tmp_pos))
    init_cov_[1] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial covariance YY");
  private_nh_.param("initial_cov_aa", tmp_pos, init_cov_[2]);
  if(!std::isnan(tmp_pos))
    init_cov_[2] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial covariance AA");
    std::cout<<"init_pose_x,y,z:"<< init_pose_[0]<<","<< init_pose_[1]<<","<< init_pose_[2]<<std::endl;
    std::cout<<"init_cov_0,1,2:"<<init_cov_[0]<<","<<init_cov_[1]<<","<<init_cov_[2]<<std::endl;
}
#pragma endregion

#pragma region LaserScanCallback
void AmclNode::laserReceived(const sensor_msgs::LaserScanConstPtr& laser_scan)
{
  std::string laser_scan_frame_id = stripSlash(laser_scan->header.frame_id);
  //用于判定是否长时间未接收到雷达数据
  last_laser_received_ts_ = ros::Time::now();
  //如果没有地图对象，将直接退出。
  if( map_ == NULL ) {
    return;
  }
  //TODO step 1 : 收到了map，开始加锁保护
  boost::recursive_mutex::scoped_lock lr(configuration_mutex_);
  int laser_index = -1;

  // TODO step 2: 检查是否有base_link -> laser 的转换，支持多个激光雷达？
  //find()函数用于查找具有给定键值k 的元素。如果找到该元素，则返回指向该元素的迭代器。否则，它返回一个指向map末尾的迭代器，即map :: end()。
  if(frame_to_laser_.find(laser_scan_frame_id) == frame_to_laser_.end())
  {
    // 如果没有该转换
    ROS_INFO("Setting up laser %d (frame_id=%s)\n", (int)frame_to_laser_.size(), laser_scan_frame_id.c_str());
    //lasers_记录下当前构建的雷达对象
    lasers_.push_back(new AMCLLaser(*laser_));
    //lasers_update_标记雷达的更新状态
    lasers_update_.push_back(true);
    laser_index = frame_to_laser_.size();

    geometry_msgs::PoseStamped ident;
    ident.header.frame_id = laser_scan_frame_id;
    ident.header.stamp = ros::Time();
    tf2::toMsg(tf2::Transform::getIdentity(), ident.pose);

    geometry_msgs::PoseStamped laser_pose;
    try
    {
      this->tf_->transform(ident, laser_pose, base_frame_id_);
    }
    catch(tf2::TransformException& e)
    {
      ROS_ERROR("Couldn't transform from %s to %s, "
                "even though the message notifier is in use",
                laser_scan_frame_id.c_str(),
                base_frame_id_.c_str());
      return;
    }

    pf_vector_t laser_pose_v;
    laser_pose_v.v[0] = laser_pose.pose.position.x;
    laser_pose_v.v[1] = laser_pose.pose.position.y;
    // laser mounting angle gets computed later -> set to 0 here!
    laser_pose_v.v[2] = 0;
    lasers_[laser_index]->SetLaserPose(laser_pose_v);
    ROS_INFO("Received laser's pose wrt robot: %.3f %.3f %.3f",laser_pose_v.v[0],laser_pose_v.v[1],laser_pose_v.v[2]);
     //通过一个string到int的map建立其雷达坐标ID到雷达对象在lasers_中的对应关系
    frame_to_laser_[laser_scan_frame_id] = laser_index;
  } 
  else 
  {
    // 如果有到laser的转换，获得laser的索引
    laser_index = frame_to_laser_[laser_scan_frame_id];
  }

  //TODO step 3: 获得机器人base在odom坐标系中的位姿（获得当前激光帧时）
  //获取接收到雷达数据时刻的里程计位姿。如果无法获取则报错退出
  pf_vector_t pose;
  if(!getOdomPose(latest_odom_pose_, pose.v[0], pose.v[1], pose.v[2], laser_scan->header.stamp, base_frame_id_))
  {
    ROS_ERROR("Couldn't determine robot's pose associated with laser scan");
    return;
  }

  pf_vector_t delta = pf_vector_zero();
  //TODO step 4: 如果不是第一帧激光扫描，看运动变化是否超过阈值，超过了则更新（第一帧是指更新了地图或者更新初始位姿）
  //*pf_init_通过前一帧和后一帧的pose信息，用于判断是否进行更新，如果不更新lasers_update_[i]=false,影响Step 5.2步和Step6
  if(pf_init_)
  {
    // 计算位姿的变化, 存储到delta 中
    //delta = pf_vector_coord_sub(pose, pf_odom_pose_);
    delta.v[0] = pose.v[0] - pf_odom_pose_.v[0];
    delta.v[1] = pose.v[1] - pf_odom_pose_.v[1];
    delta.v[2] = angle_diff(pose.v[2], pf_odom_pose_.v[2]);

    //! important: 是否需要更新滤波器
    bool update = fabs(delta.v[0]) > d_thresh_ || fabs(delta.v[1]) > d_thresh_ || fabs(delta.v[2]) > a_thresh_;
    update = update || m_force_update;
    m_force_update=false;

    // 设置laser更新的标志位
    if(update)
      for(unsigned int i=0; i < lasers_update_.size(); i++)
        lasers_update_[i] = true;
  }

  bool force_publication = false;
  //TODO step 5.1: 如果是第一帧激光扫描则初始化一些参数
  //* pf_init第一帧和默认为false，然后自身经过if函数，将其赋值为true，使其下一帧激光雷达数据到来后，可以进行Step4
  if(!pf_init_)
  {
    // 最近一次滤波器更新时的位姿，用到odom位姿
    pf_odom_pose_ = pose;

    // 滤波器现在初始化了
    pf_init_ = true;

    // 应该更新传感器数据
    for(unsigned int i=0; i < lasers_update_.size(); i++)
      lasers_update_[i] = true;

    force_publication = true;

    resample_count_ = 0;
  }
  
  //TODO step 5.2: 如果已经初始化并且已经运动了，更新运动模型
  else if(pf_init_ && lasers_update_[laser_index])
  {
    AMCLOdomData odata;
    odata.pose = pose;     //这一帧雷达对应的小车的位姿，在odom坐标系下
    odata.delta = delta;   //前一位姿与这一帧位姿的偏差。
    // 这是amcl_odom.cpp中最重要的一个函数，实现了用运动模型来更新现有的每一个粒子的位姿（这里得到的只是当前时刻的先验位姿）
    odom_->UpdateAction(pf_, (AMCLSensorData*)&odata);
  }

  bool resampled = false;
  //TODO step 6: 已经运动了，根据激光的扫描数据更新滤波器
  //*默认赋值为true，自身经过if函数，将其赋值为false。
  if(lasers_update_[laser_index])
  {
    AMCLLaserData ldata;
    ldata.sensor = lasers_[laser_index];
    ldata.range_count = laser_scan->ranges.size();

    // To account for lasers that are mounted upside-down, we determine the
    // min, max, and increment angles of the laser in the base frame.
    //
    // Construct min and max angles of laser, in the base_link frame.
    tf2::Quaternion q;
    q.setRPY(0.0, 0.0, laser_scan->angle_min);
    geometry_msgs::QuaternionStamped min_q, inc_q;
    min_q.header.stamp = laser_scan->header.stamp;
    min_q.header.frame_id = stripSlash(laser_scan->header.frame_id);
    tf2::convert(q, min_q.quaternion);

    q.setRPY(0.0, 0.0, laser_scan->angle_min + laser_scan->angle_increment);
    inc_q.header = min_q.header;
    tf2::convert(q, inc_q.quaternion);
    try
    {
      tf_->transform(min_q, min_q, base_frame_id_);
      tf_->transform(inc_q, inc_q, base_frame_id_);
    }
    catch(tf2::TransformException& e)
    {
      ROS_WARN("Unable to transform min/max laser angles into base frame: %s",
               e.what());
      return;
    }

    double angle_min = tf2::getYaw(min_q.quaternion);
    double angle_increment = tf2::getYaw(inc_q.quaternion) - angle_min;

    // wrapping angle to [-pi .. pi]
    //C 库函数 double fmod(double x, double y) 返回 x 除以 y 的余数
    angle_increment = fmod(angle_increment + 5*M_PI, 2*M_PI) - M_PI;

    ROS_DEBUG("Laser %d angles in base frame: min: %.3f inc: %.3f", laser_index, angle_min, angle_increment);

    // Apply range min/max thresholds, if the user supplied them
    if(laser_max_range_ > 0.0)
      ldata.range_max = std::min(laser_scan->range_max, (float)laser_max_range_);
    else
      ldata.range_max = laser_scan->range_max;
    double range_min;
    if(laser_min_range_ > 0.0)
      range_min = std::max(laser_scan->range_min, (float)laser_min_range_);
    else
      range_min = laser_scan->range_min;

    if(ldata.range_max <= 0.0 || range_min < 0.0) {
      ROS_ERROR("range_max or range_min from laser is negative! ignore this message.");
      return; // ignore this.
    }

    // The AMCLLaserData destructor will free this memory
    ldata.ranges = new double[ldata.range_count][2];
    ROS_ASSERT(ldata.ranges);
    for(int i=0;i<ldata.range_count;i++)
    {
      // amcl doesn't (yet) have a concept of min range.  So we'll map short
      // readings to max range.
      // 激光雷达传上来的数据只标记了最大值最小值，但是没做处理，直接将原始数据传上来，
      if(laser_scan->ranges[i] <= range_min) // 这里将最小值当最大值处理，因为在类似likelihood_field模型中，会直接将最大值丢弃
        ldata.ranges[i][0] = ldata.range_max;
      else if(laser_scan->ranges[i] > ldata.range_max)
        ldata.ranges[i][0] = std::numeric_limits<decltype(ldata.range_max)>::max();
      else
        ldata.ranges[i][0] = laser_scan->ranges[i];
      // Compute bearing
      ldata.ranges[i][1] = angle_min + (i * angle_increment);
    }

    // step 通过判断前面设置调用pf_update_sensor函数处理对应的测量模型，这里三个测量模型在《概率机器人》的第六章,
    // pf_update_sensor实现对所有粒子更新权重，并归一化、计算长期似然和短期似然,UpdateSensor接口完成粒子滤波器的测量更新。
    lasers_[laser_index]->UpdateSensor(pf_, (AMCLSensorData*)&ldata);

    lasers_update_[laser_index] = false;

    pf_odom_pose_ = pose;
    // resample_interval_次激光雷达回调之后进行粒子重采样，resample_interval_默认为2
    if(!(++resample_count_ % resample_interval_))
    {
      //按照一定的规则重采样粒子，包括前面说的失效恢复、粒子权重等，然后放入到kdtree，暂时先理解成关于位姿的二叉树，
      //然后进行聚类，得到均值和方差等信息,应该是相近的一堆粒子融合成一个粒子了，没必要维持太多相近的
      pf_update_resample(pf_);
      resampled = true;
    }

    pf_sample_set_t* set = pf_->sets + pf_->current_set;
    ROS_INFO("Num samples: %d\n", set->sample_count);

    //m_force_update默认都是false，只有通过request_nomotion_update话题信息，该标志位才为true
    if (!m_force_update)
    {
      // 将新粒子发布到全局坐标系下，一般是map
      geometry_msgs::PoseArray cloud_msg;
      cloud_msg.header.stamp = ros::Time::now();
      cloud_msg.header.frame_id = global_frame_id_;
      cloud_msg.poses.resize(set->sample_count);
      for(int i=0;i<set->sample_count;i++)
      {
        cloud_msg.poses[i].position.x = set->samples[i].pose.v[0];
        cloud_msg.poses[i].position.y = set->samples[i].pose.v[1];
        cloud_msg.poses[i].position.z = 0;
        tf2::Quaternion q;
        q.setRPY(0, 0, set->samples[i].pose.v[2]);
        tf2::convert(q, cloud_msg.poses[i].orientation);
      }
      particlecloud_pub_.publish(cloud_msg);
    }
  }

  //TODO step 7： 如果需要重采样或者强制发布
  //* resampled默认为false，只有当Step 6中if(!(++resample_count_ % resample_interval_))满足条件，resampled才为true
  //* force_publication 默认为false,当进行Step5.1时，force_publication才为true
  if(resampled || force_publication)
  {
    double max_weight = 0.0;
    int max_weight_hyp = -1;
    std::vector<amcl_hyp_t> hyps;
    hyps.resize(pf_->sets[pf_->current_set].cluster_count);
    //TODO step 7.1: 遍历所有粒子簇，找出权重均值最大的簇，其平均位姿就是我们要求的机器人后验位姿，到此一次循环已经所有完成
    for(int hyp_count = 0;hyp_count < pf_->sets[pf_->current_set].cluster_count; hyp_count++)
    {
      double weight;
      pf_vector_t pose_mean;
      pf_matrix_t pose_cov;
      if (!pf_get_cluster_stats(pf_, hyp_count, &weight, &pose_mean, &pose_cov))
      {
        ROS_ERROR("Couldn't get stats on cluster %d", hyp_count);
        break;
      }

      hyps[hyp_count].weight = weight;
      hyps[hyp_count].pf_pose_mean = pose_mean;
      hyps[hyp_count].pf_pose_cov = pose_cov;

      if(hyps[hyp_count].weight > max_weight)
      {
        max_weight = hyps[hyp_count].weight;
        max_weight_hyp = hyp_count;
      }
    }
    //TODO step 7.2： 将位姿(后验位姿xyz)，粒子集合，协方差矩阵 (和表示三个转角信息的6D协方差矩阵) 等进行更新、发布
    if(max_weight > 0.0)
    {
      ROS_DEBUG("Max weight pose: %.3f %.3f %.3f",
                 hyps[max_weight_hyp].pf_pose_mean.v[0],
                hyps[max_weight_hyp].pf_pose_mean.v[1],
                hyps[max_weight_hyp].pf_pose_mean.v[2]);

      /*
         puts("");
         pf_matrix_fprintf(hyps[max_weight_hyp].pf_pose_cov, stdout, "%6.3f");
         puts("");
       */
      //TODO step 7.2.1 用于发布的位姿和分布
      geometry_msgs::PoseWithCovarianceStamped p;
      // 填入 header
      p.header.frame_id = global_frame_id_;
      p.header.stamp = laser_scan->header.stamp;
      // 拷贝位置
      p.pose.pose.position.x = hyps[max_weight_hyp].pf_pose_mean.v[0];
      p.pose.pose.position.y = hyps[max_weight_hyp].pf_pose_mean.v[1];

      tf2::Quaternion q;
      q.setRPY(0, 0, hyps[max_weight_hyp].pf_pose_mean.v[2]);
      tf2::convert(q, p.pose.pose.orientation);
      //TODO step 7.2.2 拷贝数据到协方差， 从3D 转换成6D
      pf_sample_set_t* set = pf_->sets + pf_->current_set;
      for(int i=0; i<2; i++)
      {
        for(int j=0; j<2; j++)
        {
          // Report the overall filter covariance, rather than the
          // covariance for the highest-weight cluster
          //p.covariance[6*i+j] = hyps[max_weight_hyp].pf_pose_cov.m[i][j];
          p.pose.covariance[6*i+j] = set->cov.m[i][j];
        }
      }
      // Report the overall filter covariance, rather than the
      // covariance for the highest-weight cluster
      //p.covariance[6*5+5] = hyps[max_weight_hyp].pf_pose_cov.m[2][2];
      p.pose.covariance[6*5+5] = set->cov.m[2][2];

      /*
         printf("cov:\n");
         for(int i=0; i<6; i++)
         {
         for(int j=0; j<6; j++)
         printf("%6.3f ", p.covariance[6*i+j]);
         puts("");
         }
       */

      pose_pub_.publish(p);
      last_published_pose = p;

      ROS_DEBUG("New pose: %6.3f %6.3f %6.3f",
               hyps[max_weight_hyp].pf_pose_mean.v[0],
               hyps[max_weight_hyp].pf_pose_mean.v[1],
               hyps[max_weight_hyp].pf_pose_mean.v[2]);

      //? map->base转换中 除去odom->base得到map->odom，最后发布的是map->odom转换。
      geometry_msgs::PoseStamped odom_to_map;
      try
      {
        tf2::Quaternion q;
        q.setRPY(0, 0, hyps[max_weight_hyp].pf_pose_mean.v[2]);
        //!  tmp_tf是base_link在global map下的坐标，即base-->map
        tf2::Transform tmp_tf(q, tf2::Vector3(hyps[max_weight_hyp].pf_pose_mean.v[0], hyps[max_weight_hyp].pf_pose_mean.v[1],0.0));

        //! tmp_tf.inverse()是输入，tmp_tf_stamped.pose是输出
        //! tmp_tf_stamped是global map原点在base_link下的坐标，即map-->base
        geometry_msgs::PoseStamped tmp_tf_stamped;
        tmp_tf_stamped.header.frame_id = base_frame_id_;
        tmp_tf_stamped.header.stamp = laser_scan->header.stamp;
        tf2::toMsg(tmp_tf.inverse(), tmp_tf_stamped.pose);
        
        //! 将global map原点在base_link下的坐标变换成global map原点在odom下的坐标
        //! 即map-->odom，相当于在odom原点看map原点的位置
        //! 这里的odom_to_map并非真的odom-->map，而是反过来map-->odom
        this->tf_->transform(tmp_tf_stamped, odom_to_map, odom_frame_id_);
      }
      catch(tf2::TransformException)
      {
        ROS_DEBUG("Failed to subtract base to odom transform");
        return;
      }

      //转换odom_to_map.pose为latest_tf_
      tf2::convert(odom_to_map.pose, latest_tf_);
      latest_tf_valid_ = true;

     //* tf_broadcast_默认值为true，true才能发布map到odom的坐标变换
      if (tf_broadcast_ == true)
      {
        // We want to send a transform that is good up until a
        // tolerance time so that odom can be used
        ros::Time transform_expiration = (laser_scan->header.stamp + transform_tolerance_);
        geometry_msgs::TransformStamped tmp_tf_stamped;
        tmp_tf_stamped.header.frame_id = global_frame_id_;
        tmp_tf_stamped.header.stamp = transform_expiration;
        tmp_tf_stamped.child_frame_id = odom_frame_id_; 
        
        //! tmp_tf_stamped这个变换是odom原点在map坐标系的坐标，即odom-->map
        tf2::convert(latest_tf_.inverse(), tmp_tf_stamped.transform);

        this->tfb_->sendTransform(tmp_tf_stamped);
        sent_first_transform_ = true;
      }
    }
    else
    {
      ROS_ERROR("No pose!");
    }
  }
 
  //TODO step 8 如果上次tf有效
  //* latest_tf_valid_默认值为false，只有当粒子滤波最大权重大于0时(Step 7.2)时，latest_tf_valid_=true
  else if(latest_tf_valid_)
  {
    //TODO step 8.1
    //* tf_broadcast_默认值为true，true才能发布map到odom的坐标变换
    if (tf_broadcast_ == true)
    {
      ros::Time transform_expiration = (laser_scan->header.stamp +transform_tolerance_);
      geometry_msgs::TransformStamped tmp_tf_stamped;
      tmp_tf_stamped.header.frame_id = global_frame_id_;
      tmp_tf_stamped.header.stamp = transform_expiration;
      tmp_tf_stamped.child_frame_id = odom_frame_id_;
      tf2::convert(latest_tf_.inverse(), tmp_tf_stamped.transform);
      this->tfb_->sendTransform(tmp_tf_stamped);
    }

    //TODO step 8.2 保存最近的位姿到参数服务器
    ros::Time now = ros::Time::now();
    if((save_pose_period.toSec() > 0.0) &&  (now - save_pose_last_time) >= save_pose_period)
    {
      this->savePoseToServer();
      save_pose_last_time = now;
    }
  }

  diagnosic_updater_.update();
}

bool AmclNode::getOdomPose(geometry_msgs::PoseStamped& odom_pose,
                      double& x, double& y, double& yaw,
                      const ros::Time& t, const std::string& f)   //f参数为base_link
{
  // Get the robot's pose
  geometry_msgs::PoseStamped ident;
  ident.header.frame_id = stripSlash(f);
  ident.header.stamp = t;
  tf2::toMsg(tf2::Transform::getIdentity(), ident.pose);
  try
  {
    this->tf_->transform(ident, odom_pose, odom_frame_id_);
  }
  catch(tf2::TransformException e)
  {
    ROS_WARN("Failed to compute odom pose, skipping scan (%s)", e.what());
    return false;
  }
  x = odom_pose.pose.position.x;
  y = odom_pose.pose.position.y;
  yaw = tf2::getYaw(odom_pose.pose.orientation);
  return true;
}


void AmclNode::savePoseToServer()
{
  // We need to apply the last transform to the latest odom pose to get
  // the latest map pose to store.  We'll take the covariance from
  // last_published_pose.
  tf2::Transform odom_pose_tf2;
  tf2::convert(latest_odom_pose_.pose, odom_pose_tf2);
  tf2::Transform map_pose = latest_tf_.inverse() * odom_pose_tf2;

  double yaw = tf2::getYaw(map_pose.getRotation());

  ROS_DEBUG("Saving pose to server. x: %.3f, y: %.3f", map_pose.getOrigin().x(), map_pose.getOrigin().y() );

  private_nh_.setParam("initial_pose_x", map_pose.getOrigin().x());
  private_nh_.setParam("initial_pose_y", map_pose.getOrigin().y());
  private_nh_.setParam("initial_pose_a", yaw);
  private_nh_.setParam("initial_cov_xx",
                                  last_published_pose.pose.covariance[6*0+0]);
  private_nh_.setParam("initial_cov_yy",
                                  last_published_pose.pose.covariance[6*1+1]);
  private_nh_.setParam("initial_cov_aa",
                                  last_published_pose.pose.covariance[6*5+5]);
}
#pragma endregion

#pragma region other
void AmclNode::reconfigureCB(AMCLConfig &config, uint32_t level)
{
  boost::recursive_mutex::scoped_lock cfl(configuration_mutex_);

  //we don't want to do anything on the first call
  //which corresponds to startup
  if(first_reconfigure_call_)
  {
    first_reconfigure_call_ = false;
    default_config_ = config;
    return;
  }

  if(config.restore_defaults) {
    config = default_config_;
    //avoid looping
    config.restore_defaults = false;
  }

  d_thresh_ = config.update_min_d;
  a_thresh_ = config.update_min_a;

  resample_interval_ = config.resample_interval;

  laser_min_range_ = config.laser_min_range;
  laser_max_range_ = config.laser_max_range;

  gui_publish_period = ros::Duration(1.0/config.gui_publish_rate);
  save_pose_period = ros::Duration(1.0/config.save_pose_rate);

  transform_tolerance_.fromSec(config.transform_tolerance);

  max_beams_ = config.laser_max_beams;
  alpha1_ = config.odom_alpha1;
  alpha2_ = config.odom_alpha2;
  alpha3_ = config.odom_alpha3;
  alpha4_ = config.odom_alpha4;
  alpha5_ = config.odom_alpha5;

  z_hit_ = config.laser_z_hit;
  z_short_ = config.laser_z_short;
  z_max_ = config.laser_z_max;
  z_rand_ = config.laser_z_rand;
  sigma_hit_ = config.laser_sigma_hit;
  lambda_short_ = config.laser_lambda_short;
  laser_likelihood_max_dist_ = config.laser_likelihood_max_dist;

  if(config.laser_model_type == "beam")
    laser_model_type_ = LASER_MODEL_BEAM;
  else if(config.laser_model_type == "likelihood_field")
    laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD;
  else if(config.laser_model_type == "likelihood_field_prob")
    laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD_PROB;

  if(config.odom_model_type == "diff")
    odom_model_type_ = ODOM_MODEL_DIFF;
  else if(config.odom_model_type == "omni")
    odom_model_type_ = ODOM_MODEL_OMNI;
  else if(config.odom_model_type == "diff-corrected")
    odom_model_type_ = ODOM_MODEL_DIFF_CORRECTED;
  else if(config.odom_model_type == "omni-corrected")
    odom_model_type_ = ODOM_MODEL_OMNI_CORRECTED;

  if(config.min_particles > config.max_particles)
  {
    ROS_WARN("You've set min_particles to be greater than max particles, this isn't allowed so they'll be set to be equal.");
    config.max_particles = config.min_particles;
  }

  min_particles_ = config.min_particles;
  max_particles_ = config.max_particles;
  alpha_slow_ = config.recovery_alpha_slow;
  alpha_fast_ = config.recovery_alpha_fast;
  tf_broadcast_ = config.tf_broadcast;

  do_beamskip_= config.do_beamskip;
  beam_skip_distance_ = config.beam_skip_distance;
  beam_skip_threshold_ = config.beam_skip_threshold;

  //  清除排队的laser对象，为了更新他们的参数
  lasers_.clear();
  lasers_update_.clear();
  frame_to_laser_.clear();

  if( pf_ != NULL )
  {
    pf_free( pf_ );
    pf_ = NULL;
  }
  pf_ = pf_alloc(min_particles_, max_particles_,
                 alpha_slow_, alpha_fast_,
                 (pf_init_model_fn_t)AmclNode::uniformPoseGenerator,
                 (void *)map_);
  pf_set_selective_resampling(pf_, selective_resampling_);
  pf_err_ = config.kld_err;
  pf_z_ = config.kld_z;
  pf_->pop_err = pf_err_;
  pf_->pop_z = pf_z_;

  // 初始化滤波器
  pf_vector_t pf_init_pose_mean = pf_vector_zero();
  pf_init_pose_mean.v[0] = last_published_pose.pose.pose.position.x;
  pf_init_pose_mean.v[1] = last_published_pose.pose.pose.position.y;
  pf_init_pose_mean.v[2] = tf2::getYaw(last_published_pose.pose.pose.orientation);
  pf_matrix_t pf_init_pose_cov = pf_matrix_zero();
  pf_init_pose_cov.m[0][0] = last_published_pose.pose.covariance[6*0+0];
  pf_init_pose_cov.m[1][1] = last_published_pose.pose.covariance[6*1+1];
  pf_init_pose_cov.m[2][2] = last_published_pose.pose.covariance[6*5+5];
  pf_init(pf_, pf_init_pose_mean, pf_init_pose_cov);
  pf_init_ = false;

  // 实例化传感器对象
  // Odometry
  delete odom_;
  odom_ = new AMCLOdom();
  ROS_ASSERT(odom_);
  odom_->SetModel( odom_model_type_, alpha1_, alpha2_, alpha3_, alpha4_, alpha5_ );

  // Laser
  delete laser_;
  laser_ = new AMCLLaser(max_beams_, map_);
  ROS_ASSERT(laser_);
  if(laser_model_type_ == LASER_MODEL_BEAM)
    laser_->SetModelBeam(z_hit_, z_short_, z_max_, z_rand_,
                         sigma_hit_, lambda_short_, 0.0);
  else if(laser_model_type_ == LASER_MODEL_LIKELIHOOD_FIELD_PROB){
    ROS_INFO("Initializing likelihood field model; this can take some time on large maps...");
    laser_->SetModelLikelihoodFieldProb(z_hit_, z_rand_, sigma_hit_,
					laser_likelihood_max_dist_,
					do_beamskip_, beam_skip_distance_,
					beam_skip_threshold_, beam_skip_error_threshold_);
    ROS_INFO("Done initializing likelihood field model with probabilities.");
  }
  else if(laser_model_type_ == LASER_MODEL_LIKELIHOOD_FIELD){
    ROS_INFO("Initializing likelihood field model; this can take some time on large maps...");
    laser_->SetModelLikelihoodField(z_hit_, z_rand_, sigma_hit_,
                                    laser_likelihood_max_dist_);
    ROS_INFO("Done initializing likelihood field model.");
  }

  odom_frame_id_ = stripSlash(config.odom_frame_id);
  base_frame_id_ = stripSlash(config.base_frame_id);
  global_frame_id_ = stripSlash(config.global_frame_id);

  delete laser_scan_filter_;
  laser_scan_filter_ =
          new tf2_ros::MessageFilter<sensor_msgs::LaserScan>(*laser_scan_sub_,
                                                             *tf_,
                                                             odom_frame_id_,
                                                             100,
                                                             nh_);
  laser_scan_filter_->registerCallback(boost::bind(&AmclNode::laserReceived,
                                                   this, _1));
  // 其实调用了handleInitialPoseMessage，处理初始化位姿
  initial_pose_sub_ = nh_.subscribe("initialpose", 2, &AmclNode::initialPoseReceived, this);
}

void AmclNode::runFromBag(const std::string &in_bag_fn, bool trigger_global_localization)
{
  rosbag::Bag bag;
  bag.open(in_bag_fn, rosbag::bagmode::Read);
  std::vector<std::string> topics;
  topics.push_back(std::string("tf"));
  std::string scan_topic_name = "base_scan"; // TODO determine what topic this actually is from ROS
  topics.push_back(scan_topic_name);
  rosbag::View view(bag, rosbag::TopicQuery(topics));

  ros::Publisher laser_pub = nh_.advertise<sensor_msgs::LaserScan>(scan_topic_name, 100);
  ros::Publisher tf_pub = nh_.advertise<tf2_msgs::TFMessage>("/tf", 100);

  // Sleep for a second to let all subscribers connect
  ros::WallDuration(1.0).sleep();

  ros::WallTime start(ros::WallTime::now());

  // Wait for map
  while (ros::ok())
  {
    {
      boost::recursive_mutex::scoped_lock cfl(configuration_mutex_);
      if (map_)
      {
        ROS_INFO("Map is ready");
        break;
      }
    }
    ROS_INFO("Waiting for the map...");
    ros::getGlobalCallbackQueue()->callAvailable(ros::WallDuration(1.0));
  }

  if (trigger_global_localization)
  {
    std_srvs::Empty empty_srv;
    globalLocalizationCallback(empty_srv.request, empty_srv.response);
  }

  BOOST_FOREACH(rosbag::MessageInstance const msg, view)
  {
    if (!ros::ok())
    {
      break;
    }

    // Process any ros messages or callbacks at this point
    ros::getGlobalCallbackQueue()->callAvailable(ros::WallDuration());

    tf2_msgs::TFMessage::ConstPtr tf_msg = msg.instantiate<tf2_msgs::TFMessage>();
    if (tf_msg != NULL)
    {
      tf_pub.publish(msg);
      for (size_t ii=0; ii<tf_msg->transforms.size(); ++ii)
      {
        tf_->setTransform(tf_msg->transforms[ii], "rosbag_authority");
      }
      continue;
    }

    sensor_msgs::LaserScan::ConstPtr base_scan = msg.instantiate<sensor_msgs::LaserScan>();
    if (base_scan != NULL)
    {
      laser_pub.publish(msg);
      laser_scan_filter_->add(base_scan);
      if (bag_scan_period_ > ros::WallDuration(0))
      {
        bag_scan_period_.sleep();
      }
      continue;
    }

    ROS_WARN_STREAM("Unsupported message type" << msg.getTopic());
  }

  bag.close();

  double runtime = (ros::WallTime::now() - start).toSec();
  ROS_INFO("Bag complete, took %.1f seconds to process, shutting down", runtime);

  const geometry_msgs::Quaternion & q(last_published_pose.pose.pose.orientation);
  double yaw, pitch, roll;
  tf2::Matrix3x3(tf2::Quaternion(q.x, q.y, q.z, q.w)).getEulerYPR(yaw,pitch,roll);
  ROS_INFO("Final location %.3f, %.3f, %.3f with stamp=%f",
            last_published_pose.pose.pose.position.x,
            last_published_pose.pose.pose.position.y,
            yaw, last_published_pose.header.stamp.toSec()
            );

  ros::shutdown();
}

void AmclNode::checkLaserReceived(const ros::TimerEvent& event)
{
  ros::Duration d = ros::Time::now() - last_laser_received_ts_;
  if(d > laser_check_interval_)
  {
    ROS_WARN("No laser scan received (and thus no pose updates have been published) for %f seconds.  Verify that data is being published on the %s topic.",
             d.toSec(),
             ros::names::resolve(scan_topic_).c_str());
  }
}

AmclNode::~AmclNode()
{
  delete dsrv_;
  freeMapDependentMemory();
  delete laser_scan_filter_;
  delete laser_scan_sub_;
  // TODO: delete everything allocated in constructor
}

pf_vector_t AmclNode::uniformPoseGenerator(void* arg)
{
  map_t* map = (map_t*)arg;
#if NEW_UNIFORM_SAMPLING
  unsigned int rand_index = drand48() * free_space_indices.size();
  std::pair<int,int> free_point = free_space_indices[rand_index];
  pf_vector_t p;
  p.v[0] = MAP_WXGX(map, free_point.first);
  p.v[1] = MAP_WYGY(map, free_point.second);
  p.v[2] = drand48() * 2 * M_PI - M_PI;
#else
  double min_x, max_x, min_y, max_y;

  min_x = (map->size_x * map->scale)/2.0 - map->origin_x;
  max_x = (map->size_x * map->scale)/2.0 + map->origin_x;
  min_y = (map->size_y * map->scale)/2.0 - map->origin_y;
  max_y = (map->size_y * map->scale)/2.0 + map->origin_y;

  pf_vector_t p;

  ROS_DEBUG("Generating new uniform sample");
  for(;;)
  {
    p.v[0] = min_x + drand48() * (max_x - min_x);
    p.v[1] = min_y + drand48() * (max_y - min_y);
    p.v[2] = drand48() * 2 * M_PI - M_PI;
    // Check that it's a free cell
    int i,j;
    i = MAP_GXWX(map, p.v[0]);
    j = MAP_GYWY(map, p.v[1]);
    if(MAP_VALID(map,i,j) && (map->cells[MAP_INDEX(map,i,j)].occ_state == -1))
      break;
  }
#endif
  return p;
}

bool AmclNode::globalLocalizationCallback(std_srvs::Empty::Request& req,std_srvs::Empty::Response& res)
{
  if( map_ == NULL ) {
    return true;
  }
  boost::recursive_mutex::scoped_lock gl(configuration_mutex_);
  ROS_INFO("Initializing with uniform distribution");
  pf_init_model(pf_, (pf_init_model_fn_t)AmclNode::uniformPoseGenerator,
                (void *)map_);
  ROS_INFO("Global initialisation done!");
  pf_init_ = false;
  return true;
}

// force nomotion updates (amcl updating without requiring motion)
bool AmclNode::nomotionUpdateCallback(std_srvs::Empty::Request& req, std_srvs::Empty::Response& res)
{
	m_force_update = true;
	//ROS_INFO("Requesting no-motion update");
	return true;
}

void AmclNode::standardDeviationDiagnostics(diagnostic_updater::DiagnosticStatusWrapper& diagnostic_status)
{
  double std_x = sqrt(last_published_pose.pose.covariance[6*0+0]);
  double std_y = sqrt(last_published_pose.pose.covariance[6*1+1]);
  double std_yaw = sqrt(last_published_pose.pose.covariance[6*5+5]);

  diagnostic_status.add("std_x", std_x);
  diagnostic_status.add("std_y", std_y);
  diagnostic_status.add("std_yaw", std_yaw);
  diagnostic_status.add("std_warn_level_x", std_warn_level_x_);
  diagnostic_status.add("std_warn_level_y", std_warn_level_y_);
  diagnostic_status.add("std_warn_level_yaw", std_warn_level_yaw_);

  if (std_x > std_warn_level_x_ || std_y > std_warn_level_y_ || std_yaw > std_warn_level_yaw_)
  {
    diagnostic_status.summary(diagnostic_msgs::DiagnosticStatus::WARN, "Too large");
  }
  else
  {
    diagnostic_status.summary(diagnostic_msgs::DiagnosticStatus::OK, "OK");
  }
}
#pragma endregion