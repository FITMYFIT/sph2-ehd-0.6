#include "ControlSPH.h"

CControlSPH::CControlSPH()
{

	 _RunMod=1;//RunMod=0 Implict=0; RumMod=1 Explicit

	 _StartStep=0;//0-Start from timeStep 0=0;esle-Start from StartStep
	
	 _VSL=0;

	 _SPHAV=0;

	 _SPHPV=0;//物理粘性控制参数，0：牛顿粘性 1：幂律型粘性

	 _SPHAS=0;

	 _SPHST=0;

	 //_SPHAD=0;

	 _SPHIPCS=0;


	 _PerdBnd=0;//0 不算周期性边界，1 计算周期性边界（用于Poiseuille流等）

	 _XSPHEpsilon=0;

	 _Cs=0;//声速

	 _DeltaTCoeff=0;

	 _DeltaT=0;

	 _FinalTime=0;

	 _DensRenormSteps=0;//多少步密度重构一次

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

	 _ADDelta=0;//Artificial diffusion (人工耗散项，A Colagrossi)
	
	 _Xcrmin=0;
	 _Xcrmax=0;
	 _Ycrmin=0;
	 _Ycrmax=0;
	 _Zcrmin=0;
	 _Zcrmax=0;

	 //监测变量必须先置0
	_CSPMIflag1=0;//CSPM算法修正核函数插值时的分母，1 代表已经计算，0 还没计算
	_CSPMIflag2=0;//CSPM算法修正核函数插值时的分母，1 代表已经计算，0 还没计算

}

CControlSPH::~CControlSPH()
{
}
