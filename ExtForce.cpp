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
		ID=Region._CalList[i];//�������Ǵ�0��ʼ��
		PID=Region._PtList[ID]._PID;

		if(Region._ControlSPH._RunMod==1)// ��ʽ����
		{
			Region._PtList[ID]._du+=Region._ExtForceList[PID-1]->_fx;
			Region._PtList[ID]._dv+=Region._ExtForceList[PID-1]->_fy;
		}

		else//��ʽ����
		{
			ID2=Region._PtList[ID]._ID2;

			Region._SPHbu[ID2-1]+=Region._ExtForceList[PID-1]->_fx;
			Region._SPHbv[ID2-1]+=Region._ExtForceList[PID-1]->_fy;
		}
	}
}

void CExtForce::Solve2( CRegion &Region,double TimeSteps )//����Damp���̵�����ʩ�ӣ�XY Hu2012 Equ.(13)
{
	CSPHPt * SPHPtPtr;
	CPart * PartPtr;
	unsigned int ID,ID2,PID;
	double Time;
	double zita;//X Y Hu Equ.��13��������ʩ�ӳ˵�ϵ��
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
	
	//ʩ������
	for(unsigned int i=0;i!=Region._CalList.size();++i)
	{
		ID=Region._CalList[i];//�������Ǵ�0��ʼ��
		PID=Region._PtList[ID]._PID;

		if(Region._ControlSPH._RunMod==1)// ��ʽ����
		{
			Region._PtList[ID]._du+=zita*Region._ExtForceList[PID-1]->_fx;
			Region._PtList[ID]._dv+=zita*Region._ExtForceList[PID-1]->_fy;
		}

		else//��ʽ����
		{
			ID2=Region._PtList[ID]._ID2;

			Region._SPHbu[ID2-1]+=zita*Region._ExtForceList[PID-1]->_fx;
			Region._SPHbv[ID2-1]+=zita*Region._ExtForceList[PID-1]->_fy;
		}
	}
}
