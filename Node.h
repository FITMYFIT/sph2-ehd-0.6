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

#ifndef _NODE_H_
#define _NODE_H_

#include <vector>

enum enNODELOCALTYPE {enLeftBotm,enBotm,enRigtBotm,enRigt,enRigtUper,enUper,enLeftUper,enLeft,enIn};

class CNode
{
public:

	CNode();

	~CNode();

	enNODELOCALTYPE _LocalType;


	unsigned int _Iflag;//Iflag=0 不在计算域；Iflag=1 在计算域内

	unsigned int _ID;

	double _x;

	double _y;

	double _z;

	//the surface normal
	double _nx,_ny,_nz;
	double _Nnx,_Nny,_Nnz;	
	double _Curvature;

	double _VectMod;//法向模值

	std::vector<unsigned int> _MshIDList;//这个节点相邻的网格的编号，每个节点最多有4个相邻网格

	//double _m;	//mass

	//double _rho;//density

	//double  _rho0;

	////pressure
	//double _p;

	//double _u;//x component
	//double _v;//y component
	//double _w;//z component


	//double _h;//smoothing length
	//double _H;

	//double _r;
	//double _rr;

	//double _Volume;


	//double _uXSPHCoef;//XSPH法求解出的粒子移动速度的修正量
	//double _vXSPHCoef;
	//double _wXSPHCoef;

	//double _C0;//初始色值，定义在这重复了，待修改了计算表张程序后去掉
	//double _C;//插值后的色值

	//double _CSPMCoef;//CSPM修正函数插值时的分母

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
