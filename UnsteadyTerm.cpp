#include "UnsteadyTerm.h"

CUnsteadyTerm::CUnsteadyTerm()
{
}

CUnsteadyTerm::~CUnsteadyTerm()
{
}

void CUnsteadyTerm::Solve(CRegion &Region)
{
	CSPHPt * SPHPtPtr;
	CKnl * KnlPtr;


	unsigned int IDi,IDj;
	unsigned int ID;
	unsigned int ID2;
	double DeltaT;

	DeltaT=Region._ControlSPH._DeltaT;

	//2.2 计算源项和系数项中的非稳态项
	for(unsigned int i=0;i<Region._CalList.size();i++)
	{
		ID=Region._CalList[i];

		ID2=Region._PtList[ID]._ID2;// 实际上ID2=i+1

		//源项
		Region._SPHbu[ID2-1]+=Region._PtList[ID]._u/DeltaT;
		Region._SPHbv[ID2-1]+=Region._PtList[ID]._v/DeltaT;

		//主对角元
		Region._SPHAm[ID2-1]._Ele+=1/DeltaT;
	}
}
