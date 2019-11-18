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

	//2.2 ����Դ���ϵ�����еķ���̬��
	for(unsigned int i=0;i<Region._CalList.size();i++)
	{
		ID=Region._CalList[i];

		ID2=Region._PtList[ID]._ID2;// ʵ����ID2=i+1

		//Դ��
		Region._SPHbu[ID2-1]+=Region._PtList[ID]._u/DeltaT;
		Region._SPHbv[ID2-1]+=Region._PtList[ID]._v/DeltaT;

		//���Խ�Ԫ
		Region._SPHAm[ID2-1]._Ele+=1/DeltaT;
	}
}
