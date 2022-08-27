/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2008, 2013, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Eitan Marder-Eppstein
 *********************************************************************/
#include <costmap_2d/observation_buffer.h>

#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <tf2_sensor_msgs/tf2_sensor_msgs.h>
#include <sensor_msgs/point_cloud2_iterator.h>

using namespace std;
using namespace tf2;

namespace costmap_2d
{
ObservationBuffer::ObservationBuffer(string topic_name, double observation_keep_time, double expected_update_rate,
                                     double min_obstacle_height, double max_obstacle_height, double obstacle_range,
                                     double raytrace_range, tf2_ros::Buffer& tf2_buffer, string global_frame,
                                     string sensor_frame, double tf_tolerance) :
    tf2_buffer_(tf2_buffer), observation_keep_time_(observation_keep_time), expected_update_rate_(expected_update_rate),
    last_updated_(ros::Time::now()), global_frame_(global_frame), sensor_frame_(sensor_frame), topic_name_(topic_name),
    min_obstacle_height_(min_obstacle_height), max_obstacle_height_(max_obstacle_height),
    obstacle_range_(obstacle_range), raytrace_range_(raytrace_range), tf_tolerance_(tf_tolerance)
{
}

ObservationBuffer::~ObservationBuffer()
{
}

bool ObservationBuffer::setGlobalFrame(const std::string new_global_frame)
{
  ros::Time transform_time = ros::Time::now();
  std::string tf_error;

  geometry_msgs::TransformStamped transformStamped;
  if (!tf2_buffer_.canTransform(new_global_frame, global_frame_, transform_time, ros::Duration(tf_tolerance_), &tf_error))
  {
    ROS_ERROR("Transform between %s and %s with tolerance %.2f failed: %s.", new_global_frame.c_str(),
              global_frame_.c_str(), tf_tolerance_, tf_error.c_str());
    return false;
  }

  list<Observation>::iterator obs_it;
  for (obs_it = observation_list_.begin(); obs_it != observation_list_.end(); ++obs_it)
  {
    try
    {
      Observation& obs = *obs_it;

      geometry_msgs::PointStamped origin;
      origin.header.frame_id = global_frame_;
      origin.header.stamp = transform_time;
      origin.point = obs.origin_;

      // we need to transform the origin of the observation to the new global frame
      tf2_buffer_.transform(origin, origin, new_global_frame);
      obs.origin_ = origin.point;

      // we also need to transform the cloud of the observation to the new global frame
      tf2_buffer_.transform(*(obs.cloud_), *(obs.cloud_), new_global_frame);
    }
    catch (TransformException& ex)
    {
      ROS_ERROR("TF Error attempting to transform an observation from %s to %s: %s", global_frame_.c_str(),
                new_global_frame.c_str(), ex.what());
      return false;
    }
  }

  // now we need to update our global_frame member
  global_frame_ = new_global_frame;
  return true;
}

void ObservationBuffer::bufferCloud(const sensor_msgs::PointCloud2& cloud)
{
  geometry_msgs::PointStamped global_origin;

  // 在要填充的列表上创建新观察
  observation_list_.push_front(Observation());

  // 检查是否已明确设置起始帧，或者是否应从云中获取它
  string origin_frame = sensor_frame_ == "" ? cloud.header.frame_id : sensor_frame_;

  try
  {
    // 鉴于这些观察结果来自传感器
    geometry_msgs::PointStamped local_origin;
    local_origin.header.stamp = cloud.header.stamp;
    local_origin.header.frame_id = origin_frame;
    local_origin.point.x = 0;
    local_origin.point.y = 0;
    local_origin.point.z = 0;
    //存储传感器在global下的原始坐标点
    tf2_buffer_.transform(local_origin, global_origin, global_frame_);
    tf2::convert(global_origin.point, observation_list_.front().origin_);

    // 确保将观测缓冲区的射线扫描障碍物范围传递给观测值
    observation_list_.front().raytrace_range_ = raytrace_range_;
    observation_list_.front().obstacle_range_ = obstacle_range_;

    sensor_msgs::PointCloud2 global_frame_cloud;

    //将雷达坐标系下得到的点云数据cloud  变换到  global_frame下
    tf2_buffer_.transform(cloud, global_frame_cloud, global_frame_);
    global_frame_cloud.header.stamp = cloud.header.stamp;

    // 我们需要删除点云中低于或高于我们的高度阈值的观测值
   // std::cout<<"cloud_:"<<observation_list_.front().cloud_->data.size()<<","<<"point_step:"<< global_frame_cloud.point_step<<std::endl;
    sensor_msgs::PointCloud2& observation_cloud = *(observation_list_.front().cloud_);
    observation_cloud.height = global_frame_cloud.height;
    observation_cloud.width = global_frame_cloud.width;
    observation_cloud.fields = global_frame_cloud.fields;
    observation_cloud.is_bigendian = global_frame_cloud.is_bigendian;
    observation_cloud.point_step = global_frame_cloud.point_step;
    observation_cloud.row_step = global_frame_cloud.row_step;
    observation_cloud.is_dense = global_frame_cloud.is_dense;

    /*
    std::cout<<"width:"<<observation_cloud.width <<",height:"<<observation_cloud.height<<std::endl;
    std::cout<<"point_step:"<<observation_cloud.point_step<<",row_step:"<<observation_cloud.row_step<<std::endl;
    std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    std::cout<<"fields_size:"<<observation_cloud.fields.size()<<std::endl;
    std::cout<<"fields0: "<<(int)observation_cloud.fields[0].datatype<<" "<<observation_cloud.fields[0].name<<std::endl;
    std::cout<<"fields1: "<<(int)observation_cloud.fields[1].datatype<<" "<<observation_cloud.fields[1].name<<std::endl;
    std::cout<<"fields2: "<<(int)observation_cloud.fields[2].datatype<<" "<<observation_cloud.fields[2].name<<std::endl;
    std::cout<<"fields3: "<<(int)observation_cloud.fields[3].datatype<<" "<<observation_cloud.fields[3].name<<std::endl;
    std::cout<<"fields4: "<<(int)observation_cloud.fields[4].datatype<<" "<<observation_cloud.fields[4].name<<std::endl;
    std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    */

    unsigned int cloud_size = global_frame_cloud.height*global_frame_cloud.width;
    //允许像修改容器一样修改 sensor _ msgs: : PointCloud2
    sensor_msgs::PointCloud2Modifier modifier(observation_cloud);
    modifier.resize(cloud_size);
    unsigned int point_count = 0;

    // 复制高度范围内的点
    sensor_msgs::PointCloud2Iterator<float> iter_z(global_frame_cloud, "z");
    std::vector<unsigned char>::const_iterator iter_global = global_frame_cloud.data.begin();
    std::vector<unsigned char>::const_iterator iter_global_end = global_frame_cloud.data.end();
    std::vector<unsigned char>::iterator iter_obs = observation_cloud.data.begin();
    for (; iter_global != iter_global_end; ++iter_z, iter_global += global_frame_cloud.point_step)
    {
        //std::cout<<"iter_z:"<<*iter_z<<std::endl;
        if ((*iter_z) <= max_obstacle_height_ && (*iter_z) >= min_obstacle_height_)
        {
          /*
            std::copy(fist,last,x)
            fist [IN]: 要拷贝元素的首地址
            last [IN]:要拷贝元素的最后一个元素的下一个地址
            x [OUT] : 拷贝的目的地的首地址 
          */
          std::copy(iter_global, iter_global + global_frame_cloud.point_step, iter_obs);
          iter_obs += global_frame_cloud.point_step;
          ++point_count;
        }
    }

    // resize the cloud for the number of legal points
    modifier.resize(point_count);
    observation_cloud.header.stamp = cloud.header.stamp;
    observation_cloud.header.frame_id = global_frame_cloud.header.frame_id;
    //std::cout<<"point_count:"<<point_count<<","<<"cloud_:"<<observation_list_.front().cloud_->data.size()<<std::endl;
  }
  catch (TransformException& ex)
  {
    // if an exception occurs, we need to remove the empty observation from the list
    observation_list_.pop_front();
    ROS_ERROR("TF Exception that should never happen for sensor frame: %s, cloud frame: %s, %s", sensor_frame_.c_str(),
              cloud.header.frame_id.c_str(), ex.what());
    return;
  }

  // 如果更新成功，我们希望更新上次更新的时间
  last_updated_ = ros::Time::now();

  // 从列表中删除任何过时的观察结果
  purgeStaleObservations();
}

void ObservationBuffer::purgeStaleObservations()
{
    if (!observation_list_.empty())
    {
        list<Observation>::iterator obs_it = observation_list_.begin();
        //i如果我们没有时间保持观测……那么我们只保持一次观测
        if (observation_keep_time_ == ros::Duration(0.0))
        {
            observation_list_.erase(++obs_it, observation_list_.end());
            return;
        }

        // 否则……我们将不得不通过观测值循环，看看哪些是陈旧的
        for (obs_it = observation_list_.begin(); obs_it != observation_list_.end(); ++obs_it)
        {
            Observation& obs = *obs_it;
            // 检查观察结果是否过时……如果过时，请将其和后续内容从列表中删除
            if ((last_updated_ - obs.cloud_->header.stamp) > observation_keep_time_)
            {
              observation_list_.erase(obs_it, observation_list_.end());
              return;
            }
        }
    }
}

// returns a copy of the observations
void ObservationBuffer::getObservations(vector<Observation>& observations)
{
  // first... let's make sure that we don't have any stale observations
  purgeStaleObservations();

  // now we'll just copy the observations for the caller
  list<Observation>::iterator obs_it;
  for (obs_it = observation_list_.begin(); obs_it != observation_list_.end(); ++obs_it)
  {
    //std::cout<<"hello getObservations"<<std::endl;
    observations.push_back(*obs_it);
  }
}

bool ObservationBuffer::isCurrent() const
{
  if (expected_update_rate_ == ros::Duration(0.0))
    return true;

  bool current = (ros::Time::now() - last_updated_).toSec() <= expected_update_rate_.toSec();
  if (!current)
  {
    ROS_WARN(
        "The %s observation buffer has not been updated for %.2f seconds, and it should be updated every %.2f seconds.",
        topic_name_.c_str(), (ros::Time::now() - last_updated_).toSec(), expected_update_rate_.toSec());
  }
  return current;
}

void ObservationBuffer::resetLastUpdated()
{
  last_updated_ = ros::Time::now();
}
}  // namespace costmap_2d

