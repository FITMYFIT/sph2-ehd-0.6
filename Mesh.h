/******************************************************************************
!
! Purpose:to extend various types of particles
!
! Description:contain the basic parameters of particles
!
! Notes:
!
!******************************************************************************
!
!$Id: BASEPT,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _MESH_H_
#define _MESH_H_

#include "Node.h"
#include "BasePt.h"
#include <vector>

class CMesh
{
public:

	CMesh();

	~CMesh();

	//enLOCALTYPE _LocalType;

	unsigned int _Localx;//x方向是否在边缘1，第一层；2，第2层；3，第3层；0内部
	unsigned int _Localy;



	unsigned int _Iflag;//Iflag=0 不在计算域；Iflag=1 在计算域内

	unsigned int _ID;

	unsigned int _PID;

	unsigned int _ID2;//进入计算域的粒子编号，每步不同

	unsigned int _ID3;//进入表张计算域的粒子编号，每步不同


	double _x;

	double _y;

	double _z;



	double _h;//smoothing length
	double _H;

	double _r;
	double _rr;

	double _Volume;

	double _C0;//初始色值，定义在这重复了，待修改了计算表张程序后去掉
	double _C;//插值后的色值

	//the surface normal
	double _nx,_ny,_nz;
	double _Nnx,_Nny,_Nnz;	
	double _Curvature;

	CNode * _MeshNodList[4];//组成该网格的4个节点

	unsigned int _NebrMeshID[4];//该网格的4个方向的相邻网格的ID，按左右下上排列，如果缺哪个，则将该ID置0

	double _m;	//mass

	double _rho;//density

	double  _rho0;

	//pressure
	double _p;

	double _u;//x component
	double _v;//y component
	double _w;//z component

	//应用Remesh算法
	double _Mometu;//x方向的动量
	double _Mometv;


	double _CSPMCoef;//CSPM修正函数插值时的分母

	//double _CSPMAxx;//CSPM修正函数梯度时的系数矩阵
	//double _CSPMAxy;
	//double _CSPMAxz;
	//double _CSPMAyx;
	//double _CSPMAyy;
	//double _CSPMAyz;
	//double _CSPMAzx;
	//double _CSPMAzy;
	//double _CSPMAzz;

	//double _CSPMminLambda;//CSPM系数矩阵的最小特征值

private:
};



#endif

/******************************************************************************
!
!
!******************************************************************************/
