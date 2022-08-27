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
 *         David V. Lu!!
 *********************************************************************/
#include <algorithm>
#include <costmap_2d/inflation_layer.h>
#include <costmap_2d/costmap_math.h>
#include <costmap_2d/footprint.h>
#include <boost/thread.hpp>
#include <pluginlib/class_list_macros.h>

PLUGINLIB_EXPORT_CLASS(costmap_2d::InflationLayer, costmap_2d::Layer)

using costmap_2d::LETHAL_OBSTACLE; //254
using costmap_2d::INSCRIBED_INFLATED_OBSTACLE;   //253
using costmap_2d::NO_INFORMATION;  //255

namespace costmap_2d
{

  InflationLayer::InflationLayer()
    : resolution_(0)
    , inflation_radius_(0)
    , inscribed_radius_(0)
    , weight_(0)
    , inflate_unknown_(false)
    , cell_inflation_radius_(0)
    , cached_cell_inflation_radius_(0)
    , dsrv_(NULL)
    , seen_(NULL)
    , cached_costs_(NULL)
    , cached_distances_(NULL)
    , last_min_x_(-std::numeric_limits<float>::max())
    , last_min_y_(-std::numeric_limits<float>::max())
    , last_max_x_(std::numeric_limits<float>::max())
    , last_max_y_(std::numeric_limits<float>::max())
  {
    inflation_access_ = new boost::recursive_mutex();
  }
  
  #pragma region onInitialize
  void InflationLayer::onInitialize()
  {
    {
      boost::unique_lock < boost::recursive_mutex > lock(*inflation_access_);
      ros::NodeHandle nh("~/" + name_), g_nh;
      current_ = true;
      if (seen_)
        delete[] seen_;
      seen_ = NULL;
      seen_size_ = 0;
      need_reinflation_ = false;

      dynamic_reconfigure::Server<costmap_2d::InflationPluginConfig>::CallbackType cb = boost::bind(
          &InflationLayer::reconfigureCB, this, _1, _2);

      if (dsrv_ != NULL){
        dsrv_->clearCallback();
        dsrv_->setCallback(cb);
      }
      else
      {
        dsrv_ = new dynamic_reconfigure::Server<costmap_2d::InflationPluginConfig>(ros::NodeHandle("~/" + name_));
        dsrv_->setCallback(cb);
      }
    }
    matchSize();
  }

  void InflationLayer::reconfigureCB(costmap_2d::InflationPluginConfig &config, uint32_t level)
  {
    setInflationParameters(config.inflation_radius, config.cost_scaling_factor);

    if (enabled_ != config.enabled || inflate_unknown_ != config.inflate_unknown) {
      enabled_ = config.enabled;
      inflate_unknown_ = config.inflate_unknown;
      need_reinflation_ = true;
    }
  }

  void InflationLayer::matchSize()
  {
    boost::unique_lock < boost::recursive_mutex > lock(*inflation_access_);
    costmap_2d::Costmap2D* costmap = layered_costmap_->getCostmap();
    resolution_ = costmap->getResolution();
    //cellDistance函数功能:max(0.0,inflation_radius_/resolution_),函数在costmap_2d.cpp中
    cell_inflation_radius_ = cellDistance(inflation_radius_);
   // std::cout<<"matchSize:"<<cell_inflation_radius_<<",inscribed_radius_:"<<inscribed_radius_<<std::endl;
    computeCaches();
    unsigned int size_x = costmap->getSizeInCellsX(), size_y = costmap->getSizeInCellsY();
    
    // std::cout<<"(inflation_layer)--->"<<costmap->getSizeInCellsX()<<","<<costmap->getSizeInCellsY()<<","
    //                                                               <<costmap->getResolution()<<","<<costmap->getOriginX()<<","
    //                                                               <<costmap->getOriginY()<<std::endl;
    if (seen_)
      delete[] seen_;
    seen_size_ = size_x * size_y;
    seen_ = new bool[seen_size_];
  }

  void InflationLayer::computeCaches()
  {
    // cell_inflation_radius_为零说明不用膨胀了
    if (cell_inflation_radius_ == 0)
      return;

    // 基于膨胀半径，计算cached_distances_和cached_costs_，两者为下面膨胀计算的参照物
    if (cell_inflation_radius_ != cached_cell_inflation_radius_)
    {
      deleteKernels();
      // cached_distances_和cached_costs_的行数都是cell_inflation_radius_+2
      cached_costs_ = new unsigned char*[cell_inflation_radius_ + 2];
      cached_distances_ = new double*[cell_inflation_radius_ + 2];

      for (unsigned int i = 0; i <= cell_inflation_radius_ + 1; ++i)
      {
        // cached_distances_和cached_costs_的列数也是cell_inflation_radius_+2
        cached_costs_[i] = new unsigned char[cell_inflation_radius_ + 2];
        cached_distances_[i] = new double[cell_inflation_radius_ + 2];
        for (unsigned int j = 0; j <= cell_inflation_radius_ + 1; ++j)
        {
          //每个元素的值为它到(0,0)点的距离
          cached_distances_[i][j] = hypot(i, j);
        }
      }
      //设置cached_cell_inflation_radius_，所以第二次再进入的时候就不再执行这个if结构
      cached_cell_inflation_radius_ = cell_inflation_radius_;
    }

    for (unsigned int i = 0; i <= cell_inflation_radius_ + 1; ++i)
    {
      for (unsigned int j = 0; j <= cell_inflation_radius_ + 1; ++j)
      {
        cached_costs_[i][j] = computeCost(cached_distances_[i][j]);
        //std::cout<<"cached_costs_["<<i<<"]["<<j<<"]:"<<(float)cached_costs_[i][j]<<std::endl;
      }
    }
    // std::cout<<"cached_costs_[0][1]:"<<(float)cached_costs_[0][1]<<std::endl;
    //std::cout<<"+++++++++++++++++++++++++++++++++"<<std::endl;
  }

  void InflationLayer::deleteKernels()
  {
    if (cached_distances_ != NULL)
    {
      for (unsigned int i = 0; i <= cached_cell_inflation_radius_ + 1; ++i)
      {
        if (cached_distances_[i])
          delete[] cached_distances_[i];
      }
      if (cached_distances_)
        delete[] cached_distances_;
      cached_distances_ = NULL;
    }

    if (cached_costs_ != NULL)
    {
      for (unsigned int i = 0; i <= cached_cell_inflation_radius_ + 1; ++i)
      {
        if (cached_costs_[i])
          delete[] cached_costs_[i];
      }
      delete[] cached_costs_;
      cached_costs_ = NULL;
    }
  }
  #pragma endregion 

  #pragma region updateBounds
  // need_reinflation_默认false，更新bound这里和其他两层的主要区别是，膨胀层在传入的bound的值的基础上，通过inflation_radius_再次扩张
  void InflationLayer::updateBounds(double robot_x, double robot_y, double robot_yaw, double* min_x, double* min_y, double* max_x, double* max_y)
  {
    if (need_reinflation_)
    {
      last_min_x_ = *min_x;
      last_min_y_ = *min_y;
      last_max_x_ = *max_x;
      last_max_y_ = *max_y;
      // For some reason when I make these -<double>::max() it does not
      // work with Costmap2D::worldToMapEnforceBounds(), so I'm using
      // -<float>::max() instead.
      *min_x = -std::numeric_limits<float>::max();
      *min_y = -std::numeric_limits<float>::max();
      *max_x = std::numeric_limits<float>::max();
      *max_y = std::numeric_limits<float>::max();
      need_reinflation_ = false;
    }
    else
    {
      double tmp_min_x = last_min_x_;
      double tmp_min_y = last_min_y_;
      double tmp_max_x = last_max_x_;
      double tmp_max_y = last_max_y_;
      last_min_x_ = *min_x;
      last_min_y_ = *min_y;
      last_max_x_ = *max_x;
      last_max_y_ = *max_y;
      *min_x = std::min(tmp_min_x, *min_x) - inflation_radius_;
      *min_y = std::min(tmp_min_y, *min_y) - inflation_radius_;
      *max_x = std::max(tmp_max_x, *max_x) + inflation_radius_;
      *max_y = std::max(tmp_max_y, *max_y) + inflation_radius_;
    }
    //std::cout<<"inf-->min_x:"<<*min_x<<",max_x:"<<*max_x<<",min_y:"<<*min_y<<",max_y:"<<*max_y<<std::endl;
  }
  #pragma endregion

  #pragma region updateCosts
  // 用指针master_array指向主地图，并获取主地图的尺寸，确认seen_数组被正确设置
  void InflationLayer::updateCosts(costmap_2d::Costmap2D& master_grid, int min_i, int min_j, int max_i, int max_j)
  {
    boost::unique_lock < boost::recursive_mutex > lock(*inflation_access_);
    if (cell_inflation_radius_ == 0)
      return;

    //循环前膨胀层的list应该为空
    ROS_ASSERT_MSG(inflation_cells_.empty(), "The inflation list must be empty at the beginning of inflation");

    unsigned char* master_array = master_grid.getCharMap();
    unsigned int size_x = master_grid.getSizeInCellsX(), size_y = master_grid.getSizeInCellsY();
    
    //seen_默认为空
    if (seen_ == NULL) {
      ROS_WARN("InflationLayer::updateCosts(): seen_ array is NULL");
      seen_size_ = size_x * size_y;
      //这里分配了一个size和master map一样的bool数组
      seen_ = new bool[seen_size_];
    }
    else if (seen_size_ != size_x * size_y)
    {
      ROS_WARN("InflationLayer::updateCosts(): seen_ array size is wrong");
      delete[] seen_;
      seen_size_ = size_x * size_y;
      seen_ = new bool[seen_size_];
    }
    memset(seen_, false, size_x * size_y * sizeof(bool));

    // 边缘膨胀
    min_i -= cell_inflation_radius_;
    min_j -= cell_inflation_radius_;
    max_i += cell_inflation_radius_;
    max_j += cell_inflation_radius_;
    
    min_i = std::max(0, min_i);
    min_j = std::max(0, min_j);
    max_i = std::min(int(size_x), max_i);
    max_j = std::min(int(size_y), max_j);
    //std::cout<<"min_i:"<<min_i<<",min_j:"<<min_j<<",max_i:"<<max_i<<",max_j:"<<max_j<<std::endl;
    /*
      std::map<double, std::vector<CellData> > inflation_cells_;
      它是一个映射，由浮点数→CellData数组，CellData这个类定义在inflation_layer.h中，是专门用来记录当前cell的索引和与它最近的障碍物的索引的。
      这个映射以距离障碍物的距离为标准给bound内的cell归类，为膨胀做准备；自然距离为0对应的cell即障碍物cell本身，
      目前得到的inflation_cells_只包括障碍物本身
    */
    std::vector<CellData>& obs_bin = inflation_cells_[0.0];
    for (int j = min_j; j < max_j; j++)
    {
      for (int i = min_i; i < max_i; i++)
      {
        int index = master_grid.getIndex(i, j);
        unsigned char cost = master_array[index];
        if (cost == LETHAL_OBSTACLE)
        {//index,x_,y_,src_x_,src_y_
          obs_bin.push_back(CellData(index, i, j, i, j));
        }
      }
    }

    //接下来开始由目前的inflation_cells_进行传播，直到将所有的cell填充进去
    std::map<double, std::vector<CellData> >::iterator bin;
    //二重循环，外层遍历“键”，即距障碍物的距离；内层遍历“值”，即cell
    for (bin = inflation_cells_.begin(); bin != inflation_cells_.end(); ++bin)
    {
      //map键值对映射，键是距离最近的障碍物的距离，值是cell
      //second是“键”“值”对中的“值”，即遍历次数为celldata的个数
      for (int i = 0; i < bin->second.size(); ++i)
      {
        //在该“double键”下记录迭代到的celldata
        const CellData& cell = bin->second[i];
        //记录该celldata的索引index
        unsigned int index = cell.index_;

        // 如果已经访问过，跳过
        if (seen_[index])
        {
          continue;
        }

        seen_[index] = true;
        
        //录当前cell的坐标和离它最近的cell的坐标
        unsigned int mx = cell.x_;
        unsigned int my = cell.y_;
        unsigned int sx = cell.src_x_;
        unsigned int sy = cell.src_y_;
        
        //该cell与最近障碍物的距离来确定该cell的cost，这里称为新cost。
        unsigned char cost = costLookup(mx, my, sx, sy);
        unsigned char old_cost = master_array[index];                              
        /*
            ① 当原cost为NO_INFORMATION：
                    inflate_unknown_为true：当新cost>0，设为新cost；
                    inflate_unknown_为false：新cost≥INSCRIBED_INFLATED_OBSTACLE，设为新cost,默认为false
                区别就在于，如果inflate unknown，则当膨胀到主地图上的未知区域，只要有cost，就覆盖它；而当inflate unknown关闭，
                当膨胀到主地图上的未知区域时，只有新cost是障碍物，才覆盖它，否则维持未知。后者膨胀的范围更窄！
            ② 否则，取二者中大的cost
        */
        if (old_cost == NO_INFORMATION && (inflate_unknown_ ? (cost > FREE_SPACE) : (cost >= INSCRIBED_INFLATED_OBSTACLE)))
          master_array[index] = cost;
        else
          master_array[index] = std::max(old_cost, cost);

        // attempt to put the neighbors of the current cell onto the inflation list
        if (mx > 0)//左
          enqueue(index - 1, mx - 1, my, sx, sy);
        if (my > 0)//下
          enqueue(index - size_x, mx, my - 1, sx, sy);
        if (mx < size_x - 1)//右
          enqueue(index + 1, mx + 1, my, sx, sy);
        if (my < size_y - 1)//上
          enqueue(index + size_x, mx, my + 1, sx, sy);
      }
    }

    inflation_cells_.clear();
  }

  /**
   * @brief  Given an index of a cell in the costmap, place it into a list pending for obstacle inflation
   * @param  grid The costmap
   * @param  index The index of the cell
   * @param  mx The x coordinate of the cell (can be computed from the index, but saves time to store it)
   * @param  my The y coordinate of the cell (can be computed from the index, but saves time to store it)
   * @param  src_x The x index of the obstacle point inflation started at
   * @param  src_y The y index of the obstacle point inflation started at
   */
  inline void InflationLayer::enqueue(unsigned int index, unsigned int mx, unsigned int my, unsigned int src_x, unsigned int src_y)
  {
    if (!seen_[index])
    {
      // we compute our distance table one cell further than the inflation radius dictates so we can make the check below
      double distance = distanceLookup(mx, my, src_x, src_y);

      //找它和最近的障碍物的距离，如果超过了阈值，则直接返回
      //如果不超过阈值，就把它也加入inflation_cells_数组里
      if (distance > cell_inflation_radius_)
        return;

      // push the cell data onto the inflation list and mark
      inflation_cells_[distance].push_back(CellData(index, mx, my, src_x, src_y));
    }
  }
  #pragma endregion

  #pragma region other
  void InflationLayer::onFootprintChanged()
  {
    inscribed_radius_ = layered_costmap_->getInscribedRadius();
    cell_inflation_radius_ = cellDistance(inflation_radius_);
    computeCaches();
    need_reinflation_ = true;

    // ROS_WARN("InflationLayer::onFootprintChanged(): num footprint points: %lu,"
    //           " inscribed_radius_ = %.3f, inflation_radius_ = %.3f",
    //           layered_costmap_->getFootprint().size(), inscribed_radius_, inflation_radius_);
  }

  void InflationLayer::setInflationParameters(double inflation_radius, double cost_scaling_factor)
  {
    if (weight_ != cost_scaling_factor || inflation_radius_ != inflation_radius)
    {
      // Lock here so that reconfiguring the inflation radius doesn't cause segfaults
      // when accessing the cached arrays
      boost::unique_lock < boost::recursive_mutex > lock(*inflation_access_);

      inflation_radius_ = inflation_radius;
      cell_inflation_radius_ = cellDistance(inflation_radius_);
      //std::cout<<"setInflationParameters:"<<cell_inflation_radius_<<",inscribed_radius_:"<<inscribed_radius_<<std::endl;
      weight_ = cost_scaling_factor;
      need_reinflation_ = true;
      computeCaches();
    }
  }
  #pragma endregion


}  // namespace costmap_2d
