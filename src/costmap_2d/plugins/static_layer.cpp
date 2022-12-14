/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2008, 2013, Willow Garage, Inc.
 *  Copyright (c) 2015, Fetch Robotics, Inc.
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
 *         David V. Lu!!
 *********************************************************************/
#include <costmap_2d/static_layer.h>
#include <costmap_2d/costmap_math.h>
#include <pluginlib/class_list_macros.h>
#include <tf2/convert.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>

PLUGINLIB_EXPORT_CLASS(costmap_2d::StaticLayer, costmap_2d::Layer)

using costmap_2d::NO_INFORMATION;
using costmap_2d::LETHAL_OBSTACLE;
using costmap_2d::FREE_SPACE;

namespace costmap_2d
{

  StaticLayer::StaticLayer() : dsrv_(NULL) {}

  StaticLayer::~StaticLayer()
  {
    if (dsrv_)
      delete dsrv_;
  }

  #pragma region onInitialize
  void StaticLayer::onInitialize()
  {
    ros::NodeHandle nh("~/" + name_), g_nh;
    current_ = true;

    // global_frame_???global_costmap??????map, local_costmap??????odom_???????????????
    global_frame_ = layered_costmap_->getGlobalFrameID();

    std::string map_topic;
    nh.param("map_topic", map_topic, std::string("map"));
    nh.param("first_map_only", first_map_only_, false);
    nh.param("subscribe_to_updates", subscribe_to_updates_, false);

    nh.param("track_unknown_space", track_unknown_space_, true);
    nh.param("use_maximum", use_maximum_, false);

    int temp_lethal_threshold, temp_unknown_cost_value;
    // ??????????????????????????????100
    nh.param("lethal_cost_threshold", temp_lethal_threshold, int(100));
    // ???????????????????????????-1
    nh.param("unknown_cost_value", temp_unknown_cost_value, int(-1));
    
    nh.param("trinary_costmap", trinary_costmap_, true);

    lethal_threshold_ = std::max(std::min(temp_lethal_threshold, 100), 0);   //?????????100
    unknown_cost_value_ = temp_unknown_cost_value;  //?????????-1

    // ???static_layer.h???????????? ros::Subscriber map_sub_
    if (map_sub_.getTopic() != ros::names::resolve(map_topic))
    {
      // ???map server????????????????????????????????????
      ROS_INFO("Requesting the map...");
      map_sub_ = g_nh.subscribe(map_topic, 1, &StaticLayer::incomingMap, this);
      map_received_ = false;
      has_updated_data_ = false;

      ros::Rate r(10);
      //???incomingMap??????????????????map_received_ ??????true??? has_updated_data_ ??????true
      //??????map_received_?????????false??????????????????????????????????????????????????????????????????????????????true?????????????????????
      while (!map_received_ && g_nh.ok())
      {
        ros::spinOnce();
        r.sleep();
      }

      ROS_WARN("Received a %d X %d map at %f m/pix", getSizeInCellsX(), getSizeInCellsY(), getResolution());

      //??????subscribe_to_updates_?????????false
      if (subscribe_to_updates_)
      {
          ROS_INFO("Subscribing to updates");
          map_update_sub_ = g_nh.subscribe(map_topic + "_updates", 10, &StaticLayer::incomingUpdate, this);
      }
    }
    else
    {
      has_updated_data_ = true;
    }

    if (dsrv_)
    {
      delete dsrv_;
    }

    dsrv_ = new dynamic_reconfigure::Server<costmap_2d::GenericPluginConfig>(nh);
    dynamic_reconfigure::Server<costmap_2d::GenericPluginConfig>::CallbackType cb = boost::bind(
        &StaticLayer::reconfigureCB, this, _1, _2);
    dsrv_->setCallback(cb);
  }
  
  //?????????????????????
  void StaticLayer::incomingMap(const nav_msgs::OccupancyGridConstPtr& new_map)
  {
    // step 1 ??????????????????????????????????????????
    unsigned int size_x = new_map->info.width, size_y = new_map->info.height;
    ROS_INFO("Received a %d X %d map at %f m/pix", size_x, size_y, new_map->info.resolution);

    // step 2 ??????master costmap?????????????????????????????????????????????
    Costmap2D* master = layered_costmap_->getCostmap();
    #pragma region print
    // std::cout<<"**********Master Map Information:**********"<<std::endl;
    // std::cout<<"width:"<<master->getSizeInCellsX()<<std::endl;
    // std::cout<<"hight:"<<master->getSizeInCellsY()<<std::endl;
    // std::cout<<"resolution:"<<master->getResolution()<<std::endl;
    // std::cout<<"originX:"<<master->getOriginX()<<std::endl;
    // std::cout<<"originY:"<<master->getOriginY()<<std::endl;
    // std::cout<<"*********************************************"<<std::endl;
    
    // std::cout<<"**********New Map Information:**************"<<std::endl;
    // std::cout<<"width:"<<size_x<<std::endl;
    // std::cout<<"hight:"<<size_y<<std::endl;
    // std::cout<<"resolution:"<<new_map->info.resolution<<std::endl;
    // std::cout<<"originX:"<<new_map->info.origin.position.x<<std::endl;
    // std::cout<<"originY:"<<new_map->info.origin.position.y<<std::endl;
    // std::cout<<"*********************************************"<<std::endl;
    #pragma endregion
  //???layered_cosmap->isRolling==false???
    if (!layered_costmap_->isRolling() &&
        (master->getSizeInCellsX() != size_x ||
        master->getSizeInCellsY() != size_y ||
        master->getResolution() != new_map->info.resolution ||
        master->getOriginX() != new_map->info.origin.position.x ||
        master->getOriginY() != new_map->info.origin.position.y))
    {
      // ??????master costmap ???????????????????????????layered_costmap_ ??????
      ROS_INFO("Master resizing costmap to %d X %d at %f m/pix", size_x, size_y, new_map->info.resolution);
      layered_costmap_->resizeMap(size_x, size_y, new_map->info.resolution, new_map->info.origin.position.x,
                                  new_map->info.origin.position.y,
                                  true /* set size_locked to true, prevents reconfigureCb from overriding map size*/);
    }
    //???layered_cosmap->isRolling==true???
    //?????????size_x_???0,??????????????????????????????static costmap ???
    else if (size_x_ != size_x || size_y_ != size_y ||
            resolution_ != new_map->info.resolution ||
            origin_x_ != new_map->info.origin.position.x ||
            origin_y_ != new_map->info.origin.position.y)
    {
      // ?????????static costmap ?????????
      ROS_INFO("Only resizing static layer to %d X %d at %f m/pix", size_x, size_y, new_map->info.resolution);
      //??????resizeMap???costmap_2d.cpp????????????
      resizeMap(size_x, size_y, new_map->info.resolution,new_map->info.origin.position.x, new_map->info.origin.position.y);
    }

    // size_x_???size_y_??????Costmap2D?????????????????????
    // ??????size_x_???size_y_???resolution_???origin_x_???origin_y_?????????????????????map????????????
    unsigned int index = 0;

    //???????????????map????????????????????????static map???????????????costmap
    for (unsigned int i = 0; i < size_y; ++i)
    {
      for (unsigned int j = 0; j < size_x; ++j)
      {
        // ?????????????????????-1(???255),0 ???100
        unsigned char value = new_map->data[index];
        //costmap_?????????unsigned char *,????????? costmap_2d.h???
        costmap_[index] = interpretValue(value);
        ++index;
      }
    }
    
    map_frame_ = new_map->header.frame_id;  //?????????????????????
    
    // ?????????????????????????????????
    x_ = y_ = 0;
    width_ = size_x_; // ???????????????????????????(???????????????)????????????????????????int??????
    height_ = size_y_;
    map_received_ = true;
    has_updated_data_ = true;

    // ??????first_map_only_?????????????????????????????????false?????????map??????
    if (first_map_only_)
    {
      ROS_INFO("Shutting down the map subscriber. first_map_only flag is on");
      map_sub_.shutdown();
    }
  }

  unsigned char StaticLayer::interpretValue(unsigned char value)
  {
    //track_unknown_space_?????????true,unknown_cost_value_?????????-1
    if (track_unknown_space_ && value == unknown_cost_value_)
      // map?????????-1 ??????????????????
      return NO_INFORMATION;    //NO_INFORMATION==255

    else if (!track_unknown_space_ && value == unknown_cost_value_)
      return FREE_SPACE;     //FREE_SPACE=0

    // map?????????100 ??????>=100 ??????????????????????????????????????????
    else if (value >= lethal_threshold_)
      return LETHAL_OBSTACLE;    //LETHAL_OBSTACLE=254

    else if (trinary_costmap_)  //??????trinary_costmap_ ?????????true?????????????????????????????????????????????-1-->255, 0, 254-->100
      return FREE_SPACE;   

    // map?????????0-100,???????????????????????????????????????
    double scale = (double) value / lethal_threshold_;
    return scale * LETHAL_OBSTACLE;
  }

  void StaticLayer::incomingUpdate(const map_msgs::OccupancyGridUpdateConstPtr& update)
  {
    unsigned int di = 0;
    for (unsigned int y = 0; y < update->height ; y++)
    {
      unsigned int index_base = (update->y + y) * size_x_;
      for (unsigned int x = 0; x < update->width ; x++)
      {
        unsigned int index = index_base + x + update->x;
        costmap_[index] = interpretValue(update->data[di++]);
      }
    }
    x_ = update->x;
    y_ = update->y;
    width_ = update->width;
    height_ = update->height;
    has_updated_data_ = true;
  }

  void StaticLayer::reconfigureCB(costmap_2d::GenericPluginConfig &config, uint32_t level)
  {
    if (config.enabled != enabled_)
    {
      enabled_ = config.enabled;
      has_updated_data_ = true;
      x_ = y_ = 0;
      width_ = size_x_;
      height_ = size_y_;
      ROS_WARN("width:%f\n",width_);
    }
  }
  #pragma endregion

  #pragma region updateBounds
  //????????????????????????
  void StaticLayer::updateBounds(double robot_x, double robot_y, double robot_yaw, double* min_x, double* min_y,
                                double* max_x, double* max_y)
  {
     // ?????? min_x=min_y=1e30???max_x=max_y=-1e30. ????????????robot_x???robot_y,robot_yaw
    if( !layered_costmap_->isRolling() )
    {
      //has_extra_bounds_????????????false,???costmap_layer.h???????????????
      if (!map_received_ || !(has_updated_data_ || has_extra_bounds_))
        return;
    }

    // ???????????????,useExtraBounds?????????costmap_layer.cpp
    useExtraBounds(min_x, min_y, max_x, max_y);
    //std::cout<<"static2:"<<"min_x:"<<*min_x<<",min_y:"<<*min_y<<",max_x:"<<*max_x<<",max_y:"<<*max_y<<std::endl;
    
    double wx, wy;

    // ????????????costmap??????????????????????????????,mapToWorld?????????costmap_2d.cpp???
    //x_=y_=0
    mapToWorld(x_, y_, wx, wy);
    *min_x = std::min(wx, *min_x);
    *min_y = std::min(wy, *min_y);

    // ????????????costmap??????????????????????????????
    mapToWorld(x_ + width_, y_ + height_, wx, wy);
    *max_x = std::max(wx, *max_x);
    *max_y = std::max(wy, *max_y);
    // std::cout<<"*********static layer:***************"<<std::endl;
    std::cout<<"sta-->min_x:"<<*min_x<<",max_x:"<<*max_x<<",min_y:"<<*min_y<<",max_y:"<<*max_y<<std::endl;
    has_updated_data_ = false;
  }
  #pragma endregion

  #pragma region updateCosts
  //????????????????????????
  void StaticLayer::updateCosts(costmap_2d::Costmap2D& master_grid, int min_i, int min_j, int max_i, int max_j)
  {
      if (!map_received_)
        return;

      if (!layered_costmap_->isRolling())
      {
          // ??????costmap ??????????????????????????????static costmap ???cost ??????????????????master costmap
          //  static costmap ???cost ?????????incomingMap??????????????????????????????????????????costmap_2d ?????????
          //use_maximum_?????????false
          if (!use_maximum_)
            updateWithTrueOverwrite(master_grid, min_i, min_j, max_i, max_j);  //?????????costmap_layer.cpp???
          else
            updateWithMax(master_grid, min_i, min_j, max_i, max_j);  //?????????costmap_layer.cpp???
      }
      else
      {
          // ??????rolling window, master_grid??????????????????????????????????????????
          unsigned int mx, my;
          double wx, wy;
          // ???????????????????????????
          geometry_msgs::TransformStamped transform;
          try
          {
            transform = tf_->lookupTransform(map_frame_, global_frame_, ros::Time(0));
          }
          catch (tf2::TransformException ex)
          {
            ROS_ERROR("%s", ex.what());
            return;
          }
          // Copy map data given proper transformations
          tf2::Transform tf2_transform;
          tf2::convert(transform.transform, tf2_transform);
          for (unsigned int i = min_i; i < max_i; ++i)
          {
              for (unsigned int j = min_j; j < max_j; ++j)
              {
                  // Convert master_grid coordinates (i,j) into global_frame_(wx,wy) coordinates
                  layered_costmap_->getCostmap()->mapToWorld(i, j, wx, wy);
                  // Transform from global_frame_ to map_frame_
                  tf2::Vector3 p(wx, wy, 0);
                  p = tf2_transform*p;
                  // Set master_grid with cell from map
                  if (worldToMap(p.x(), p.y(), mx, my))
                  {
                    if (!use_maximum_)
                      master_grid.setCost(i, j, getCost(mx, my));
                    else
                      master_grid.setCost(i, j, std::max(getCost(mx, my), master_grid.getCost(i, j)));
                  }
              }
          }
      }
    }
  #pragma endregion

  #pragma region other
  void StaticLayer::matchSize()
  {
    // If we are using rolling costmap, the static map size is
    //   unrelated to the size of the layered costmap
    if (!layered_costmap_->isRolling())
    {
      Costmap2D* master = layered_costmap_->getCostmap();
      std::cout<<"(static_layer)--->"<<master->getSizeInCellsX()<<","<<master->getSizeInCellsY()<<","
                                                                  <<master->getResolution()<<","<<master->getOriginX()<<","
                                                                  <<master->getOriginY()<<std::endl;
      resizeMap(master->getSizeInCellsX(), master->getSizeInCellsY(), master->getResolution(),
                master->getOriginX(), master->getOriginY());
    }
  }

  void StaticLayer::activate()
  {
    onInitialize();
  }

  void StaticLayer::deactivate()
  {
    map_sub_.shutdown();
    if (subscribe_to_updates_)
      map_update_sub_.shutdown();
  }

  void StaticLayer::reset()
  {
    if (first_map_only_)
    {
      has_updated_data_ = true;
    }
    else
    {
      onInitialize();
    }
  }
  #pragma endregion

}  // namespace costmap_2d
