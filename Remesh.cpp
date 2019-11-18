#include "Remesh.h"


CRemesh::CRemesh()
{

}

CRemesh::~CRemesh()
{

}

void CRemesh::Solve( CRegion & Region )
{
	CBasePt * PtiPtr;
	CMesh * MshjPtr;
	CKnl * KnlPtr;

	//找到网格-粒子对
	_GetMshNeigb.GetMshNbl(Region);
	//_GetMshKnlList.GetMshKnlList(Region);

	//将网格上的待插值变量置0
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		MshjPtr=&Region._MeshList[i];

		MshjPtr->_m=0.0;
		MshjPtr->_rho=0.0;
		MshjPtr->_u=0.0;
		MshjPtr->_v=0.0;

		MshjPtr->_Mometu=0.0;
		MshjPtr->_Mometv=0.0;

		MshjPtr->_CSPMCoef=0.0;
	}

	//将粒子信息插值到网格上
	for (unsigned int i=0;i!=Region._PtMshPairList.size();++i)
	{
		PtiPtr=Region._PtMshPairList[i]._PtiPtr;
		MshjPtr=Region._PtMshPairList[i]._MshjPtr;
	}

	////将粒子的变量插值到背景网格上
	//for (unsigned int i=0;i!=Region._PtMshPairList.size();++i)
	//{
	//	PtiPtr=Region._PtMshPairList[i]._PtiPtr;
	//	MshjPtr=Region._PtMshPairList[i]._MshjPtr;

	//	KnlPtr=&Region._PtMshKnlList[i];

	//	MshjPtr->_rho+=PtiPtr->_m*KnlPtr->_W;

	//	MshjPtr->_u+=PtiPtr->_m/PtiPtr->_rho*PtiPtr->_u*KnlPtr->_W;
	//	MshjPtr->_v+=PtiPtr->_m/PtiPtr->_rho*PtiPtr->_v*KnlPtr->_W;

	//	MshjPtr->_CSPMCoef+=PtiPtr->_m/PtiPtr->_rho*KnlPtr->_W;
	//}

	////CSPM修正
	//for (unsigned int i=0;i!=Region._MeshList.size();++i)
	//{
	//	MshjPtr=&Region._MeshList[i];

	//	if (MshjPtr->_CSPMCoef!=0)
	//	{		
	//		MshjPtr->_rho/=MshjPtr->_CSPMCoef;

	//		MshjPtr->_u/=MshjPtr->_CSPMCoef;
	//		MshjPtr->_v/=MshjPtr->_CSPMCoef;

	//		MshjPtr->_m=MshjPtr->_rho*MshjPtr->_Volume;
	//	}
	//}

	//将网格信息转移为粒子
	std::vector<CBasePt> TempPtList;

	//将Null粒子转存
	for (unsigned int i=0;i!=Region._PtList.size();++i)
	{
		if (Region._PtList[i]._Type==enNULLPt)
		{
			TempPtList.push_back(Region._PtList[i]);
		}
	}

	Region._PtList.clear();
	unsigned int icount=0;
	for (unsigned int i=0;i!=TempPtList.size();++i)
	{
		icount++;
		TempPtList[i]._ID=icount;
		Region._PtList.push_back(TempPtList[i]);
	}
	TempPtList.clear();

	//
	CBasePt BasePt;

	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		BasePt._Type=enSPHPt;

		BasePt._x=Region._MeshList[i]._x;
		BasePt._y=Region._MeshList[i]._y;
		BasePt._m=Region._MeshList[i]._m;
		BasePt._rho=Region._MeshList[i]._rho;

		BasePt._u=Region._MeshList[i]._u;
		BasePt._v=Region._MeshList[i]._v;
		BasePt._Volume=Region._MeshList[i]._Volume;
	
		icount++;
		BasePt._ID=icount;
		BasePt._PID=2;
		
		BasePt._Cs=Region._ControlSPH._Cs;

		BasePt._h=Region._PartList[BasePt._PID-1]._HdivDp*sqrt(BasePt._Volume);
		BasePt._r=2*BasePt._h;

		BasePt._H=Region._PartList[BasePt._PID-1]._STHvsh*BasePt._h;
		BasePt._rr=2*BasePt._H;

		BasePt._rho0=1000;

		Region._PtList.push_back(BasePt);
	}


}
