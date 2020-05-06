/******************************************************************************
!
! Purpose:to control all the variables of the processes
!
! Description:contain all of the control variables of the SPH particles
!
! Notes:
!
!******************************************************************************
!
!$Id: CONTROLSPH,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _CONTROLSPH_H_
#define _CONTROLSPH_H_

#include <string>
#include <iostream>
#include <sstream>
#include <iosfwd>

class CControlSPH
{
public:

	CControlSPH();

	~CControlSPH();

	unsigned int _RunMod;//RunMod=0 Implict; RumMod=1 Explicit

	unsigned int _StartStep;//0-Start from timeStep 0;esle-Start from StartStep

	unsigned int _VSL;

	unsigned int _SPHAV;

	unsigned int _SPHPV;//物理粘性控制参数，0：牛顿粘性 1：幂律型粘性

	unsigned int _SPHAS;

	unsigned int _SPHST;

	//unsigned int _SPHAD;

	unsigned int _SPHIPCS;


	unsigned int _PerdBnd;//0 不算周期性边界，1 计算周期性边界（用于Poiseuille流等）

	double _XSPHEpsilon;//XSPH的参数

	double _Cs;//声速

	double _DeltaTCoeff;

	double _DeltaT;

	double _FinalTime;

	unsigned int _DensRenormSteps;//多少步密度重构一次

	unsigned int _OutputSteps;

	unsigned int _CellNumx;
	unsigned int _CellNumy;
	unsigned int _CellNumz;

	double _AVAlpha;
	double _AVBeta;
	double _AVEta;

	double _ASEpsilon1;
	double _ASEpsilon2;
	double _ASDeltaD;

	double _ADDelta;//Artificial diffusion (人工耗散项，A Colagrossi)

	double _Xcrmin;
	double _Xcrmax;
	double _Ycrmin;
	double _Ycrmax;
	double _Zcrmin;
	double _Zcrmax;


	//表面张力搜索域的范围（用于射流断裂计算），只在Y方向
	double _STYmin;
	double _STYmax;


	//监测变量，监测程序运行中算了哪些量

	unsigned int _CSPMIflag1;//CSPM算法修正核函数插值时的分母，1 代表已经计算，0 还没计算
	unsigned int _CSPMIflag2;//CSPM算法修正函数导数插值时的系数矩阵

	//插值用的网格信息
	double _MshMinX;
	double _MshMaxX;

	double _MshMinY;
	double _MshMaxY;

	double _MshSize;

	std::string _InfileName;//读入文件的名称


	double _PerdBndMinX;//周期性边界条件的最小x，用于Poiseuille流
	double _PerdBndMaxX;

	unsigned int _IFlagRemesh;//是否进行了Remeh

	double _TDamp;//Damp时间，XYHu 2012Equ（13）

	double _PtShftBeta;//particle shift coefficient, beta

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
