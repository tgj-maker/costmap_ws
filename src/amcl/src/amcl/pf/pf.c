/*
 *  Player - One Hell of a Robot Server
 *  Copyright (C) 2000  Brian Gerkey   &  Kasper Stoy
 *                      gerkey@usc.edu    kaspers@robotics.usc.edu
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
/**************************************************************************
 * Desc: Simple particle filter for localization.
 * Author: Andrew Howard
 * Date: 10 Dec 2002
 * CVS: $Id: pf.c 6345 2008-04-17 01:36:39Z gerkey $
 *************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "amcl/pf/pf.h"
#include "amcl/pf/pf_pdf.h"
#include "amcl/pf/pf_kdtree.h"
#include "portable_utils.hpp"
#include <float.h>

// Compute the required number of samples, given that there are k bins
// with samples in them.
static int pf_resample_limit(pf_t *pf, int k);

void copy_set(pf_sample_set_t* set_a, pf_sample_set_t* set_b);

#pragma region pf init
// Create a new filter
pf_t *pf_alloc(int min_samples, int max_samples,
               double alpha_slow, double alpha_fast,
               pf_init_model_fn_t random_pose_fn, void *random_pose_data)
{
  int i, j;
  pf_t *pf;
  pf_sample_set_t *set;   //多个位姿粒子群
  pf_sample_t *sample;  //单个粒子的位姿和权重

  srand48(time(NULL));

  /*
  说明：calloc 在内存的动态存储区中分配num个长度为size的连续空间，函数返回一个指向分配起始地址的指针；
               如果分配不成功，返回NULL。
  用法：malloc(单位：字节)：malloc(10 * sizeof(int));或malloc(40)
               calloc:calloc(10 , sizeof(int))
  注意：malloc的使用效率较高，因为calloc在返回在堆区申请的那块动态内存的起始地址之前，会将每个字节都初始化为0
  */
  pf = calloc(1, sizeof(pf_t));

  pf->min_samples = min_samples;
  pf->max_samples = max_samples;
  pf->pop_err = 0.01; // [err]是真实分布和估计分布之间的最大误差。
  pf->pop_z = 3;  // [z]是(1 -p)的上标准正态分位数，其中p是估计分布上的误差小于[err]的概率。
  // 叶节点的数量永远不会大于最大的样本数量
  pf->limit_cache = calloc(max_samples, sizeof(int));
  pf->current_set = 0;

  pf->random_pose_fn = random_pose_fn;
  pf->random_pose_data = random_pose_data;

  for (j = 0; j < 2; j++)
  {
    set = pf->sets + j;  //set类似迭代器指针，通过j进行指针位置的移动，指针类型为pf_sample_set_t

    set->sample_count = max_samples;
    //set->samples类似与指针数组，该指针刚开始指向数组的第一个元素，其中数组个数为max_samples,元素类型为pf_sample_t*
    set->samples = calloc(max_samples, sizeof(pf_sample_t)); 

    for (i = 0; i < set->sample_count; i++)
    {
      sample = set->samples + i;  //sample类似迭代器指针，通过i进行指针位置的移动，指针类型为pf_sample_t
      sample->pose.v[0] = 0.0;
      sample->pose.v[1] = 0.0;
      sample->pose.v[2] = 0.0;
      sample->weight = 1.0 / max_samples;
    }

    // 构建了三倍于样本集合尺寸的直方图对象kdtree
    set->kdtree = pf_kdtree_alloc(3 * max_samples);
    
    //初始化当前粒子群的数量为0
    set->cluster_count = 0;
    //初始化一个粒子算一类
    set->cluster_max_count = max_samples;
    //set->clusters类似与指针数组，该指针刚开始指向数组的第一个元素，其中数组个数为set->cluster_max_count,元素类型为pf_cluster_t*
    set->clusters = calloc(set->cluster_max_count, sizeof(pf_cluster_t));

    set->mean = pf_vector_zero();
    set->cov = pf_matrix_zero();
  }

  pf->w_slow = 0.0;
  pf->w_fast = 0.0;

  pf->alpha_slow = alpha_slow;
  pf->alpha_fast = alpha_fast;
  pf->dist_threshold = 0.5;
  //设定滤波器未收敛
  pf_init_converged(pf);

  return pf;
}

// Initialize the filter using a gaussian
void pf_init(pf_t *pf, pf_vector_t mean, pf_matrix_t cov)
{
  int i;
  pf_sample_set_t *set;        //多个位姿粒子群
  pf_sample_t *sample;        //单个粒子的位姿和权重
  pf_pdf_gaussian_t *pdf;   

  //通过索引current_set获取当前激活的粒子集合
  set = pf->sets + pf->current_set;

  // 创建用于自适应采样的kd树
  pf_kdtree_clear(set->kdtree);

  set->sample_count = pf->max_samples;
  
  //根据输入的均值和方差构建一个高斯分布
  pdf = pf_pdf_gaussian_alloc(mean, cov);

  // 计算新的样本姿势
  for (i = 0; i < set->sample_count; i++)
  {
    sample = set->samples + i;   //sample类似迭代器指针，通过i进行指针位置的移动，指针类型为pf_sample_t
    sample->weight = 1.0 / pf->max_samples;
    sample->pose = pf_pdf_gaussian_sample(pdf);

    // Add sample to histogram
    pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
  }

  pf->w_slow = pf->w_fast = 0.0;
  
  //释放掉临时申请的高斯分布对象pdf
  pf_pdf_gaussian_free(pdf);

  //重新计算集群统计
  pf_cluster_stats(pf, set);

  //标记为未收敛
  pf_init_converged(pf);

  return;
}

// Initialize the filter using some model
void pf_init_model(pf_t *pf, pf_init_model_fn_t init_fn, void *init_data)
{
  int i;
  pf_sample_set_t *set;
  pf_sample_t *sample;

  set = pf->sets + pf->current_set;

  // Create the kd tree for adaptive sampling
  pf_kdtree_clear(set->kdtree);

  set->sample_count = pf->max_samples;

  // Compute the new sample poses
  for (i = 0; i < set->sample_count; i++)
  {
    sample = set->samples + i;
    sample->weight = 1.0 / pf->max_samples;
    sample->pose = (*init_fn) (init_data);

    // Add sample to histogram
    pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
  }

  pf->w_slow = pf->w_fast = 0.0;

  // Re-compute cluster statistics
  pf_cluster_stats(pf, set);

  //set converged to 0
  pf_init_converged(pf);

  return;
}

void pf_init_converged(pf_t *pf){
  pf_sample_set_t *set;
  set = pf->sets + pf->current_set;
  set->converged = 0;
  pf->converged = 0;
}

// Re-compute the cluster statistics for a sample set
void pf_cluster_stats(pf_t *pf, pf_sample_set_t *set)
{
  int i, j, k, cidx;
  pf_sample_t *sample;
  pf_cluster_t *cluster;

  // Workspace
  double m[4], c[2][2];
  size_t count;
  double weight;

  //在这里就给这棵kd树打标签了，其实针对的元素是粒子sample
  pf_kdtree_cluster(set->kdtree);

  // Initialize cluster stats
  set->cluster_count = 0;

  for (i = 0; i < set->cluster_max_count; i++)
  {
    cluster = set->clusters + i;
    cluster->count = 0;
    cluster->weight = 0;
    cluster->mean = pf_vector_zero();
    cluster->cov = pf_matrix_zero();

    for (j = 0; j < 4; j++)
      cluster->m[j] = 0.0;
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        cluster->c[j][k] = 0.0;
  }
 
  // 初始化滤波器的统计特性
  count = 0;
  weight = 0.0;
  set->mean = pf_vector_zero();
  set->cov = pf_matrix_zero();
  for (j = 0; j < 4; j++)
    m[j] = 0.0;
  for (j = 0; j < 2; j++)
    for (k = 0; k < 2; k++)
      c[j][k] = 0.0;

//计算粒子群的统计特性，可以得到粒子群的均值和线性/角度方向的协方差，
//记住它是由粒子sample的权重和位姿计算得到。
  for (i = 0; i < set->sample_count; i++)
  {
    sample = set->samples + i;

    //获取粒子群的cluster标签
    cidx = pf_kdtree_get_cluster(set->kdtree, sample->pose);
    assert(cidx >= 0);
    if (cidx >= set->cluster_max_count)
      continue;
    if (cidx + 1 > set->cluster_count)
      set->cluster_count = cidx + 1;

    cluster = set->clusters + cidx;

    cluster->count += 1;
    cluster->weight += sample->weight;

    count += 1;
    weight += sample->weight;

    //计算粒子群的均值
    cluster->m[0] += sample->weight * sample->pose.v[0];
    cluster->m[1] += sample->weight * sample->pose.v[1];
    cluster->m[2] += sample->weight * cos(sample->pose.v[2]);
    cluster->m[3] += sample->weight * sin(sample->pose.v[2]);

    m[0] += sample->weight * sample->pose.v[0];
    m[1] += sample->weight * sample->pose.v[1];
    m[2] += sample->weight * cos(sample->pose.v[2]);
    m[3] += sample->weight * sin(sample->pose.v[2]);

    // 计算粒子群在线性方向上的协方差
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
      {
        cluster->c[j][k] += sample->weight * sample->pose.v[j] * sample->pose.v[k];
        c[j][k] += sample->weight * sample->pose.v[j] * sample->pose.v[k];
      }
  }

  //归一化
  for (i = 0; i < set->cluster_count; i++)
  {
    cluster = set->clusters + i;

    cluster->mean.v[0] = cluster->m[0] / cluster->weight;
    cluster->mean.v[1] = cluster->m[1] / cluster->weight;
    cluster->mean.v[2] = atan2(cluster->m[3], cluster->m[2]);

    cluster->cov = pf_matrix_zero();

    // 计算粒子群在线性方向上的协方差
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        cluster->cov.m[j][k] = cluster->c[j][k] / cluster->weight - cluster->mean.v[j] * cluster->mean.v[k];

    // Covariance in angular components; I think this is the correct
    // formula for circular statistics.
    //计算粒子群在角度方向上的协方差
    cluster->cov.m[2][2] = -2 * log(sqrt(cluster->m[2] * cluster->m[2] + cluster->m[3] * cluster->m[3]));
  }

  assert(fabs(weight) >= DBL_EPSILON);
  if (fabs(weight) < DBL_EPSILON)
  {
    printf("ERROR : divide-by-zero exception : weight is zero\n");
    return;
  }

  //计算粒子集set的统计特性：均值
  set->mean.v[0] = m[0] / weight;
  set->mean.v[1] = m[1] / weight;
  set->mean.v[2] = atan2(m[3], m[2]);


  //计算粒子集set的统计特性：协方差
  for (j = 0; j < 2; j++)
    for (k = 0; k < 2; k++)
      set->cov.m[j][k] = c[j][k] / weight - set->mean.v[j] * set->mean.v[k];

  // Covariance in angular components; I think this is the correct
  // formula for circular statistics.
  set->cov.m[2][2] = -2 * log(sqrt(m[2] * m[2] + m[3] * m[3]));

  return;
}

// Free an existing filter
void pf_free(pf_t *pf)
{
  int i;

  free(pf->limit_cache);

  for (i = 0; i < 2; i++)
  {
    free(pf->sets[i].clusters);
    pf_kdtree_free(pf->sets[i].kdtree);
    free(pf->sets[i].samples);
  }
  free(pf);

  return;
}
#pragma endregion

#pragma region pf Resample
// 用新的观测数据更新粒子滤波器，完成的是短期和长期样本均值的估计， 也就是参数w_slow和w_fast的计算
//参数pf是滤波器对象，
//参数sensor_fn则是一个函数指针实现了传感器的模型，用于计算各个粒子的概率权重。
//参数sensor_data则是用于更新的传感器数据。
void pf_update_sensor(pf_t *pf, pf_sensor_model_fn_t sensor_fn, void *sensor_data)
{
  int i;
  pf_sample_set_t *set;
  pf_sample_t *sample;
  double total;
  
  //获取当前激活的样本集合，set类似迭代器指针，通过pf->current_set进行指针位置的移动，指针类型为pf_sample_set_t
  set = pf->sets + pf->current_set;

  // 计算各个样本的概率权重， 这里sensor_fn没具体定义，具体实现会在调用pf_update_sensor时传参进来
  total = (*sensor_fn) (sensor_data, set);

  set->n_effective = 0;

  if (total > 0.0)
  {
    //归一化权重
    double w_avg=0.0;
    for (i = 0; i < set->sample_count; i++)
    {
      sample = set->samples + i;   //sample类似迭代器指针，通过i进行指针位置的移动，指针类型为pf_sample_t
      w_avg += sample->weight;
      sample->weight /= total;
      set->n_effective += sample->weight*sample->weight;
    }
    // Update running averages of likelihood of samples (Prob Rob p258)
    w_avg /= set->sample_count;
    if(pf->w_slow == 0.0)
      pf->w_slow = w_avg;
    else
      pf->w_slow += pf->alpha_slow * (w_avg - pf->w_slow);
    if(pf->w_fast == 0.0)
      pf->w_fast = w_avg;
    else
      pf->w_fast += pf->alpha_fast * (w_avg - pf->w_fast);
    //printf("w_avg: %e slow: %e fast: %e\n",
           //w_avg, pf->w_slow, pf->w_fast);
  }
  else
  {
    // Handle zero total
    for (i = 0; i < set->sample_count; i++)
    {
      sample = set->samples + i;
      sample->weight = 1.0 / set->sample_count;
    }
  }

  set->n_effective = 1.0/set->n_effective;
  return;
}

// 具体实现了重采样操作，相比于传统的粒子滤波器的重采样操作，它多了一个根据统计量w_slow,w_fast插入随机样本的过程
//参数pf用于指示进行重采样更新的滤波器对象
void pf_update_resample(pf_t *pf)
{
  int i;
  double total;
  //set_a获取当前激活的粒子群集合，set_b待切换的集合
  pf_sample_set_t *set_a, *set_b;
  pf_sample_t *sample_a, *sample_b;

  double* c;

  double w_diff;

  set_a = pf->sets + pf->current_set;
  set_b = pf->sets + (pf->current_set + 1) % 2;

  if (pf->selective_resampling != 0)
  {
    if (set_a->n_effective > 0.5*(set_a->sample_count))
    {
      // copy set a to b
      copy_set(set_a,set_b);

      // Re-compute cluster statistics
      pf_cluster_stats(pf, set_b);

      // Use the newly created sample set
      pf->current_set = (pf->current_set + 1) % 2;
      return;
    }
  }

  //建立累积概率表进行重采样
  // TODO: Replace this with a more efficient procedure
  // (e.g., http://www.network-theory.co.uk/docs/gslref/GeneralDiscreteDistributions.html)
  c = (double*)malloc(sizeof(double)*(set_a->sample_count+1));
  c[0] = 0.0;
  for(i=0;i<set_a->sample_count;i++)
    c[i+1] = c[i]+set_a->samples[i].weight;

  // 创建kd树进行自适应采样
  pf_kdtree_clear(set_b->kdtree);

  // Draw samples from set a to create set b.
  total = 0;
  set_b->sample_count = 0;

  w_diff = 1.0 - pf->w_fast / pf->w_slow;
  if(w_diff < 0.0)
    w_diff = 0.0;
  //printf("w_diff: %9.6f\n", w_diff);

 //我们重置更新粒子集合set_b，它将保存重采样后的粒子
//当粒子sample集set_b的sample数目小于粒子滤波器的sample的最大数目
  while(set_b->sample_count < pf->max_samples)
  {
    ///逐步从粒子sample集的set_b中取出sample对象
    sample_b = set_b->samples + set_b->sample_count++;

    if(drand48() < w_diff)
      //使用随机生成位姿的模型生成sample_b的位姿；
      sample_b->pose = (pf->random_pose_fn)(pf->random_pose_data);
    else  // 如果判据不生效，那么我们就根据原粒子集合所描述的样本分布进行采样
    {
      //不能(轻易地)将低方差采样器与KLD自适应采样相结合，所以我们将采取更传统的方法。
      
      double r;
      r = drand48(); //生成随机树
      for(i=0;i<set_a->sample_count;i++)
      {
        if((c[i] <= r) && (r < c[i+1]))
          break;
      }
      assert(i<set_a->sample_count);

      sample_a = set_a->samples + i;

      assert(sample_a->weight > 0);

      //把sample_a的位姿赋给sample_b
      sample_b->pose = sample_a->pose;
    }


   //重采样后的粒子赋予相同的权重，并更新统计直方图
    sample_b->weight = 1.0;
    total += sample_b->weight;
   //把粒子sample_b添加到kdtree里，其实也就是添加到粒子sample集set_b里
    pf_kdtree_insert(set_b->kdtree, sample_b->pose, sample_b->weight);

    //检查看是否有足够的粒子sample了,当满足样本限制的时候，就退出循环，结束重采样。
    //这里使用pf_resample_limit函数，输入为粒子sample集set_b的kdtree的叶子节点个数，这也代表着直方图的k个bin。
    if (set_b->sample_count > pf_resample_limit(pf, set_b->kdtree->leaf_count))
      break;
  }

  //重置平均值，以避免陷入完全的随机性。
  if(w_diff > 0.0)
    pf->w_slow = pf->w_fast = 0.0;

  // 归一化样本权重
  for (i = 0; i < set_b->sample_count; i++)
  {
    sample_b = set_b->samples + i;
    sample_b->weight /= total;
  }

  // 更新粒子群
  pf_cluster_stats(pf, set_b);

  // 交换当前粒子集
  pf->current_set = (pf->current_set + 1) % 2;

  //检查粒子群是否收敛
  pf_update_converged(pf);

  free(c);
  return;
}

// copy set a to set b
void copy_set(pf_sample_set_t* set_a, pf_sample_set_t* set_b)
{
  int i;
  double total;
  pf_sample_t *sample_a, *sample_b;

  // Clean set b's kdtree
  pf_kdtree_clear(set_b->kdtree);

  // Copy samples from set a to create set b
  total = 0;
  set_b->sample_count = 0;

  for(i = 0; i < set_a->sample_count; i++)
  {
    sample_b = set_b->samples + set_b->sample_count++;

    sample_a = set_a->samples + i;

    assert(sample_a->weight > 0);

    // Copy sample a to sample b
    sample_b->pose = sample_a->pose;
    sample_b->weight = sample_a->weight;

    total += sample_b->weight;

    // Add sample to histogram
    pf_kdtree_insert(set_b->kdtree, sample_b->pose, sample_b->weight);
  }

  // Normalize weights
  for (i = 0; i < set_b->sample_count; i++)
  {
    sample_b = set_b->samples + i;
    sample_b->weight /= total;
  }

  set_b->converged = set_a->converged;
}

int pf_update_converged(pf_t *pf)
{
  int i;
  pf_sample_set_t *set;
  pf_sample_t *sample;
  double total;

  set = pf->sets + pf->current_set;
  double mean_x = 0, mean_y = 0;

  for (i = 0; i < set->sample_count; i++){
    sample = set->samples + i;

    mean_x += sample->pose.v[0];
    mean_y += sample->pose.v[1];
  }
  mean_x /= set->sample_count;
  mean_y /= set->sample_count;

  for (i = 0; i < set->sample_count; i++){
    sample = set->samples + i;
    if(fabs(sample->pose.v[0] - mean_x) > pf->dist_threshold ||
       fabs(sample->pose.v[1] - mean_y) > pf->dist_threshold){
      set->converged = 0;
      pf->converged = 0;
      return 0;
    }
  }
  set->converged = 1;
  pf->converged = 1;
  return 1;
}

//根据直方图统计量k， 和误差边界参数ε和σ计算样本集合数量。 结合刚刚的重采样过程就是KLD_Sampling_MCL算法了。
int pf_resample_limit(pf_t *pf, int k)
{
  double a, b, c, x;
  int n;

  //在k超出预期范围的情况下返回max_samples，这不应该发生，但添加它是为了防止任何运行时错误
  if (k < 1 || k > pf->max_samples)
      return pf->max_samples;

  // 如果缓存有效，返回值，这意味着值是非零正的
  if (pf->limit_cache[k-1] > 0)
    return pf->limit_cache[k-1];

  if (k == 1)
  {
    pf->limit_cache[k-1] = pf->max_samples;
    return pf->max_samples;
  }

  a = 1;
  b = 2 / (9 * ((double) k - 1));
  c = sqrt(2 / (9 * ((double) k - 1))) * pf->pop_z;
  x = a - b + c;

  n = (int) ceil((k - 1) / (2 * pf->pop_err) * x * x * x);

  if (n < pf->min_samples)
  {
    pf->limit_cache[k-1] = pf->min_samples;
    return pf->min_samples;
  }
  if (n > pf->max_samples)
  {
    pf->limit_cache[k-1] = pf->max_samples;
    return pf->max_samples;
  }

  pf->limit_cache[k-1] = n;
  return n;
}
#pragma endregion s

#pragma region other
//设置选择性重采样的开关，amcl_node.cpp调用
void pf_set_selective_resampling(pf_t *pf, int selective_resampling)
{
  pf->selective_resampling = selective_resampling;
}

// 计算cep状态(mean 和variance) ，pf_draw.c调用
void pf_get_cep_stats(pf_t *pf, pf_vector_t *mean, double *var)
{
  int i;
  double mn, mx, my, mrr;
  pf_sample_set_t *set;
  pf_sample_t *sample;

  set = pf->sets + pf->current_set;

  mn = 0.0;
  mx = 0.0;
  my = 0.0;
  mrr = 0.0;

  for (i = 0; i < set->sample_count; i++)
  {
    sample = set->samples + i;

    mn += sample->weight;
    mx += sample->weight * sample->pose.v[0];
    my += sample->weight * sample->pose.v[1];
    mrr += sample->weight * sample->pose.v[0] * sample->pose.v[0];
    mrr += sample->weight * sample->pose.v[1] * sample->pose.v[1];
  }

  assert(fabs(mn) >= DBL_EPSILON);
  if (fabs(mn) < DBL_EPSILON)
  {
    printf("ERROR : divide-by-zero exception : mn is zero\n");
    return;
  }

  mean->v[0] = mx / mn;
  mean->v[1] = my / mn;
  mean->v[2] = 0.0;

  *var = mrr / mn - (mx * mx / (mn * mn) + my * my / (mn * mn));

  return;
}

// 获得粒子群的状态，amcl_node.cpp调用
int pf_get_cluster_stats(pf_t *pf, int clabel, double *weight,
                         pf_vector_t *mean, pf_matrix_t *cov)
{
  pf_sample_set_t *set;
  pf_cluster_t *cluster;

  set = pf->sets + pf->current_set;

  if (clabel >= set->cluster_count)
    return 0;
  cluster = set->clusters + clabel;

  *weight = cluster->weight;
  *mean = cluster->mean;
  *cov = cluster->cov;

  return 1;
}

// Update the filter with some new action 没有用到？
void pf_update_action(pf_t *pf, pf_action_model_fn_t action_fn, void *action_data)
{
  pf_sample_set_t *set;

  set = pf->sets + pf->current_set;

  (*action_fn) (action_data, set);

  return;
}
#pragma endregion 

