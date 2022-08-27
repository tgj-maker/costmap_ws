#include<costmap_2d/costmap_layer.h>

namespace costmap_2d
{

//一个bound和一个坐标，若坐标不在bound范围内，它扩张bound，使其包含坐标
void CostmapLayer::touch(double x, double y, double* min_x, double* min_y, double* max_x, double* max_y)
{
    *min_x = std::min(x, *min_x);
    *min_y = std::min(y, *min_y);
    *max_x = std::max(x, *max_x);
    *max_y = std::max(y, *max_y);
}

//用主地图的尺寸来设置该层地图的尺寸
void CostmapLayer::matchSize()
{
    Costmap2D* master = layered_costmap_->getCostmap();
    std::cout<<"(costmap_layer)--->"<<master->getSizeInCellsX()<<","<<master->getSizeInCellsY()<<","
                                                                       <<master->getResolution()<<","<<master->getOriginX()<<","
                                                                      <<master->getOriginY()<<std::endl;
    resizeMap(master->getSizeInCellsX(), master->getSizeInCellsY(), master->getResolution(),
            master->getOriginX(), master->getOriginY());
}

void CostmapLayer::clearArea(int start_x, int start_y, int end_x, int end_y, bool invert_area)
{
  unsigned char* grid = getCharMap();
  for(int x=0; x<(int)getSizeInCellsX(); x++){
    bool xrange = x>start_x && x<end_x;

    for(int y=0; y<(int)getSizeInCellsY(); y++){
      if((xrange && y>start_y && y<end_y)!=invert_area)
        continue;
      int index = getIndex(x,y);
      // 将栅格更新为未知信息
      if(grid[index]!=NO_INFORMATION){
        grid[index] = NO_INFORMATION;
      }
    }
  }
}

//将传入的bound与数据成员的值比较，如果传入的bound范围更大，则更新数据成员的值
void CostmapLayer::addExtraBounds(double mx0, double my0, double mx1, double my1)
{
    extra_min_x_ = std::min(mx0, extra_min_x_); // 1e6
    extra_max_x_ = std::max(mx1, extra_max_x_); // -1e6
    extra_min_y_ = std::min(my0, extra_min_y_);
    extra_max_y_ = std::max(my1, extra_max_y_);
    has_extra_bounds_ = true;
}

void CostmapLayer::useExtraBounds(double* min_x, double* min_y, double* max_x, double* max_y)
{
    if (!has_extra_bounds_)
        return;
    *min_x = std::min(extra_min_x_, *min_x);
    *min_y = std::min(extra_min_y_, *min_y);
    *max_x = std::max(extra_max_x_, *max_x);
    *max_y = std::max(extra_max_y_, *max_y);
    extra_min_x_ = 1e6;
    extra_min_y_ = 1e6;
    extra_max_x_ = -1e6;
    extra_max_y_ = -1e6;
    has_extra_bounds_ = false;
}

//用当前子地图数据（不包括未知cell）更新主地图对应区域，若对应cell的cost值比主地图大或主地图该cell为未知时，
//用子地图数据覆盖，否则保留主地图原数据
void CostmapLayer::updateWithMax(costmap_2d::Costmap2D& master_grid, int min_i, int min_j, int max_i, int max_j)
{
    if (!enabled_)
      return;

    unsigned char* master_array = master_grid.getCharMap();
    unsigned int span = master_grid.getSizeInCellsX();

    for (int j = min_j; j < max_j; j++)
    {
      unsigned int it = j * span + min_i;
      for (int i = min_i; i < max_i; i++)
      {
        //此时的costmap_是unsigned char *类型,代表各个层的代价值
        if (costmap_[it] == NO_INFORMATION)
        {
            it++;
            continue;
        }

        unsigned char old_cost = master_array[it];
        if (old_cost == NO_INFORMATION || old_cost < costmap_[it])
              master_array[it] = costmap_[it];
        it++;
        
      }
    }
}

//用当前子地图数据（包括未知cell）覆盖主地图对应区域
void CostmapLayer::updateWithTrueOverwrite(costmap_2d::Costmap2D& master_grid, int min_i, int min_j,
                                           int max_i, int max_j)
{
  if (!enabled_)
    return;
  unsigned char* master = master_grid.getCharMap();
  unsigned int span = master_grid.getSizeInCellsX();
  for (int j = min_j; j < max_j; j++)
  {
    unsigned int it = span*j+min_i;
    for (int i = min_i; i < max_i; i++)
    {
     //此时的costmap_是unsigned char *类型,代表各个层的代价值
      master[it] = costmap_[it];
      it++;
    }
  }
}
//用当前子地图数据（不包括未知cell）覆盖主地图对应区域
void CostmapLayer::updateWithOverwrite(costmap_2d::Costmap2D& master_grid, int min_i, int min_j, int max_i, int max_j)
{
  if (!enabled_)
    return;
  unsigned char* master = master_grid.getCharMap();
  unsigned int span = master_grid.getSizeInCellsX();

  for (int j = min_j; j < max_j; j++)
  {
    unsigned int it = span*j+min_i;
    for (int i = min_i; i < max_i; i++)
    {
      if (costmap_[it] != NO_INFORMATION)
        master[it] = costmap_[it];
      it++;
    }
  }
}

//函数用当前子地图数据（不包括未知cell）更新主地图对应区域，若主地图该cell为未知，用子地图数据覆盖；
//否则，在主地图原数据基础上+子地图数据（将进行限制避免cost值溢出）
void CostmapLayer::updateWithAddition(costmap_2d::Costmap2D& master_grid, int min_i, int min_j, int max_i, int max_j)
{
  if (!enabled_)
    return;
  unsigned char* master_array = master_grid.getCharMap();
  unsigned int span = master_grid.getSizeInCellsX();

  for (int j = min_j; j < max_j; j++)
  {
    unsigned int it = j * span + min_i;
    for (int i = min_i; i < max_i; i++)
    {
      if (costmap_[it] == NO_INFORMATION){
        it++;
        continue;
      }

      unsigned char old_cost = master_array[it];
      if (old_cost == NO_INFORMATION)
        master_array[it] = costmap_[it];
      else
      {
        int sum = old_cost + costmap_[it];
        if (sum >= costmap_2d::INSCRIBED_INFLATED_OBSTACLE)
            master_array[it] = costmap_2d::INSCRIBED_INFLATED_OBSTACLE - 1;
        else
            master_array[it] = sum;
      }
      it++;
    }
  }
}
}  // namespace costmap_2d
