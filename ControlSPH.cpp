#include "ControlSPH.h"

CControlSPH::CControlSPH()
{

	 _RunMod=1;//RunMod=0 Implict=0; RumMod=1 Explicit

	 _StartStep=0;//0-Start from timeStep 0=0;esle-Start from StartStep
	
	 _VSL=0;

	 _SPHAV=0;

	 _SPHPV=0;//����ճ�Կ��Ʋ�����0��ţ��ճ�� 1��������ճ��

	 _SPHAS=0;

	 _SPHST=0;

	 //_SPHAD=0;

	 _SPHIPCS=0;


	 _PerdBnd=0;//0 ���������Ա߽磬1 ���������Ա߽磨����Poiseuille���ȣ�

	 _XSPHEpsilon=0;

	 _Cs=0;//����

	 _DeltaTCoeff=0;

	 _DeltaT=0;

	 _FinalTime=0;

	 _DensRenormSteps=0;//���ٲ��ܶ��ع�һ��

	 _OutputSteps=0;

	 _CellNumx=0;
	 _CellNumy=0;
	 _CellNumz=0;

	 _AVAlpha=0;
	 _AVBeta=0;
	 _AVEta=0;

	 _ASEpsilon1=0;
	 _ASEpsilon2=0;
	 _ASDeltaD=0;

	 _ADDelta=0;//Artificial diffusion (�˹���ɢ�A Colagrossi)
	
	 _Xcrmin=0;
	 _Xcrmax=0;
	 _Ycrmin=0;
	 _Ycrmax=0;
	 _Zcrmin=0;
	 _Zcrmax=0;

	 //��������������0
	_CSPMIflag1=0;//CSPM�㷨�����˺�����ֵʱ�ķ�ĸ��1 �����Ѿ����㣬0 ��û����
	_CSPMIflag2=0;//CSPM�㷨�����˺�����ֵʱ�ķ�ĸ��1 �����Ѿ����㣬0 ��û����

}

CControlSPH::~CControlSPH()
{
}
