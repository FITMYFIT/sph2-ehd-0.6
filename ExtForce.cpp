#include "ExtForce.h"


CExtForce::CExtForce()
{
}

CExtForce::~CExtForce()
{
}

void CExtForce::Solve(CRegion &Region)
{
	CSPHPt * SPHPtPtr;
	CPart * PartPtr;
	unsigned int ID,ID2,PID;

	for(unsigned int i=0;i!=Region._CalList.size();++i)
	{
		ID=Region._CalList[i];//这个编号是从0开始的
		PID=Region._PtList[ID]._PID;

		if(Region._ControlSPH._RunMod==1)// 显式计算
		{
			Region._PtList[ID]._du+=Region._ExtForceList[PID-1]->_fx;
			Region._PtList[ID]._dv+=Region._ExtForceList[PID-1]->_fy;
		}

		else//隐式计算
		{
			ID2=Region._PtList[ID]._ID2;

			Region._SPHbu[ID2-1]+=Region._ExtForceList[PID-1]->_fx;
			Region._SPHbv[ID2-1]+=Region._ExtForceList[PID-1]->_fy;
		}
	}
}

void CExtForce::Solve2( CRegion &Region,double TimeSteps )//加入Damp过程的外力施加，XY Hu2012 Equ.(13)
{
	CSPHPt * SPHPtPtr;
	CPart * PartPtr;
	unsigned int ID,ID2,PID;
	double Time;
	double zita;//X Y Hu Equ.（13），外力施加乘的系数
	double TDamp;

	TDamp=Region._ControlSPH._TDamp;

	if (TDamp!=0)
	{
		Time=Region._ControlSPH._DeltaT*TimeSteps;

		if (Time<=TDamp)
		{
			zita=0.5*(sin(3.14159*(-0.5+Time/TDamp))+1);
		}

		else
		{
			zita=1;
		}
	}

	else
	{
		zita=1;
	}
	
	//施加外力
	for(unsigned int i=0;i!=Region._CalList.size();++i)
	{
		ID=Region._CalList[i];//这个编号是从0开始的
		PID=Region._PtList[ID]._PID;

		if(Region._ControlSPH._RunMod==1)// 显式计算
		{
			Region._PtList[ID]._du+=zita*Region._ExtForceList[PID-1]->_fx;
			Region._PtList[ID]._dv+=zita*Region._ExtForceList[PID-1]->_fy;
		}

		else//隐式计算
		{
			ID2=Region._PtList[ID]._ID2;

			Region._SPHbu[ID2-1]+=zita*Region._ExtForceList[PID-1]->_fx;
			Region._SPHbv[ID2-1]+=zita*Region._ExtForceList[PID-1]->_fy;
		}
	}
}
