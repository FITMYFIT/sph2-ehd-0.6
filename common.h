/***************************************************************************
 *   define some common values may be used for the whole project           *
 *   some are copied from OpenFVM                                          *
 *   2019.09.23 Richard LIU @ Beihang University                             *

 *               *
 ***************************************************************************/

#ifndef _COMMON_H_
#define _COMMON_H_

#include <math.h>
#include <float.h>

#define LOGICAL_TRUE   1
#define LOGICAL_FALSE  0
#define LOGICAL_ERROR -1

#define PI 3.14159265358979323846264
#define e 2.71828182

#define LMAX(A,B) ((A)>(B) ? (A):(B))
#define LMIN(A,B) ((A)<(B) ? (A):(B))

#define LABS(X)   ((X) < 0 ? -(X):(X))

#define LSGN(X)   ((X) < 0 ? -(1):(1))

#define ISZERO(X) (fabs(X) < 10.0 * DBL_MIN)
#define ISONE(X)  (fabs(X - 1.0) < 10.0 * DBL_EPSILON)

#define SMALL  1E-15
#define VSMALL 1E-300

#define GREAT  1E+15
#define VGREAT 1E+300

#define MAXFACES 6


//marocs for sph storage, former, not used 2019

#define NBPARTICLE 2000           //SPH模型粒子数
#define NBGHOSTPARTICLE 2000      //
#define NBPARTICLEPAIR 250000     //SPH模型粒子对数
#define NBPART 1                  //SPH模型Part数
#define NBSECTION 2               //
#define NBEXTFORCE 1              //SPH模型受到的外力数

#define NBCELLX 100 //x方向的网格数
#define NBCELLY 40  //y方向的网格数
#define MINX 0      //模型x坐标的最小值
#define MAXX 2.5e-3 //模型x坐标的最大值
#define MINY 0     //模型y坐标的最小值
#define MAXY 1.0E-3 //模型y坐标的最大值

#define MshSizeX 2.5E-5 //x方向的网格尺寸
#define MshSizeY 2.5E-5 //y方向的网格尺寸

#define MaxPtPerMsh 1   //粒子搜索时每个盒子最多放多少个粒子

#define SORalphaP 0.5   //SIMPLE算法中的压力亚松弛因子
#define SORalphaU 0.7   //SIMPLE算法中的速度U亚松弛因子
#define SORalphaV 0.7   //SIMPLE算法中的速度V亚松弛因子


#endif
