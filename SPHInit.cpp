#include "SPHInit.h"


CSPHInit::CSPHInit()
{
}

CSPHInit::~CSPHInit()
{
}

void CSPHInit::Solve(CRegion &Region,unsigned int TimeSteps)
{
  unsigned int jj;
  CSPHPt * SPHPtPtr;
  CPart * PartPtr;
  CBasePt * BasePtPtr;

  //监控参数置0
  Region._ControlSPH._CSPMIflag1=0;//监控是否计算了CSPM修正函数的分母项
  Region._ControlSPH._CSPMIflag2=0;//监控是否计算了CSPM修正函数梯度的系数项


  if(Region._ControlSPH._RunMod==1)//显式算法
	{

      for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
          BasePtPtr=&Region._PtList[i];


          BasePtPtr->_du=0.0;
          BasePtPtr->_dv=0.0;
          BasePtPtr->_dw=0.0;
          BasePtPtr->_de=0.0;
          BasePtPtr->_dh=0.0;
          BasePtPtr->_drho=0.0;

          BasePtPtr->_CSPMCoef=0.0;
          BasePtPtr->_CSPMAxx=0.0;
          BasePtPtr->_CSPMAyx=0.0;
          BasePtPtr->_CSPMAzx=0.0;
          BasePtPtr->_CSPMAxy=0.0;
          BasePtPtr->_CSPMAyy=0.0;
          BasePtPtr->_CSPMAzy=0.0;
          BasePtPtr->_CSPMAxz=0.0;
          BasePtPtr->_CSPMAyz=0.0;
          BasePtPtr->_CSPMAzz=0.0;

          BasePtPtr->_CSPMminLambda=0.0;

          BasePtPtr->_sxx=0.0;
          BasePtPtr->_sxy=0.0;
          BasePtPtr->_sxz=0.0;
          BasePtPtr->_syy=0.0;
          BasePtPtr->_syz=0.0;
          BasePtPtr->_szz=0.0;

          BasePtPtr->_AccBndx=0.0;//边界粒子对流体粒子产生的加速度
          BasePtPtr->_AccBndy=0.0;//边界粒子对流体粒子产生的加速度

          //BasePtPtr->_NumNegbor=0;//每个粒子支持域内的粒子数置0
          BasePtPtr->_r0=0.0;//截止距离，用于particle shift
          BasePtPtr->_PtShftCoefu=0.0;
          BasePtPtr->_PtShftCoefv=0.0;
          BasePtPtr->_PtShftCoefw=0.0;

          //parameters in EHD model 2019.09.12
          BasePtPtr->_eEpsilon=0.0;
          BasePtPtr->_eKappa=0.0;
          BasePtPtr->_deRho=0.0;
          BasePtPtr->_geEpsilonx=0.0;
          BasePtPtr->_geEpsilony=0.0;
          BasePtPtr->_Fex=0.0;//temp, for debug 2019.10.09
          BasePtPtr->_Fey=0.0;//temp, for debug 2019.10.09

          BasePtPtr->_mrho=BasePtPtr->_m/BasePtPtr->_rho;//2019.09.13

          if (BasePtPtr->_Type==enSPHPt||BasePtPtr->_Type==enNULLPt)
			{
              BasePtPtr->_nx=0.0;
              BasePtPtr->_ny=0.0;
              BasePtPtr->_nz=0.0;
              BasePtPtr->_Nnx=0.0;
              BasePtPtr->_Nny=0.0;
              BasePtPtr->_Nnz=0.0;
              BasePtPtr->_Curvature=0.0;
              BasePtPtr->_midcurve=0.0;
              BasePtPtr->_midcorrect=0.0;
              BasePtPtr->_C=0.0;
			}
		}
	}

  //隐式算法
  if(Region._ControlSPH._RunMod==0)
	{
      //对计算中用到的变量置初值
      for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
          BasePtPtr=&Region._PtList[i];


          BasePtPtr->_du=0.0;
          BasePtPtr->_dv=0.0;
          BasePtPtr->_dw=0.0;
          BasePtPtr->_de=0.0;
          BasePtPtr->_dh=0.0;
          BasePtPtr->_drho=0.0;

          BasePtPtr->_CSPMCoef=0.0;
          BasePtPtr->_CSPMAxx=0.0;
          BasePtPtr->_CSPMAyx=0.0;
          BasePtPtr->_CSPMAzx=0.0;
          BasePtPtr->_CSPMAxy=0.0;
          BasePtPtr->_CSPMAyy=0.0;
          BasePtPtr->_CSPMAzy=0.0;
          BasePtPtr->_CSPMAxz=0.0;
          BasePtPtr->_CSPMAyz=0.0;
          BasePtPtr->_CSPMAzz=0.0;

          BasePtPtr->_CSPMminLambda=0.0;

          BasePtPtr->_sxx=0.0;
          BasePtPtr->_sxy=0.0;
          BasePtPtr->_sxz=0.0;
          BasePtPtr->_syy=0.0;
          BasePtPtr->_syz=0.0;
          BasePtPtr->_szz=0.0;

          BasePtPtr->_AccBndx=0.0;//边界粒子对流体粒子产生的加速度
          BasePtPtr->_AccBndy=0.0;//边界粒子对流体粒子产生的加速度

          // BasePtPtr->_NumNegbor=0;//每个粒子支持域内的粒子数置0
          BasePtPtr->_r0=0.0;//截止距离，用于particle shift
          BasePtPtr->_PtShftCoefu=0.0;
          BasePtPtr->_PtShftCoefv=0.0;
          BasePtPtr->_PtShftCoefw=0.0;

          if (BasePtPtr->_Type==enSPHPt||BasePtPtr->_Type==enNULLPt)
			{
              BasePtPtr->_nx=0.0;
              BasePtPtr->_ny=0.0;
              BasePtPtr->_nz=0.0;
              BasePtPtr->_Nnx=0.0;
              BasePtPtr->_Nny=0.0;
              BasePtPtr->_Nnz=0.0;
              BasePtPtr->_Curvature=0.0;
              BasePtPtr->_midcurve=0.0;
              BasePtPtr->_midcorrect=0.0;
              BasePtPtr->_C=0.0;
			}
		}

      Region._SPHAm.resize(Region._StatDataList._InvolvedFluidNum);//只分配了主对角元，其它的需要push_back
      Region._SPHu.resize(Region._StatDataList._InvolvedFluidNum);
      Region._SPHv.resize(Region._StatDataList._InvolvedFluidNum);
      Region._SPHbu.resize(Region._StatDataList._InvolvedFluidNum);
      Region._SPHbv.resize(Region._StatDataList._InvolvedFluidNum);

      //将容器中的元素置合适的初值
      unsigned int icount=0;
      for(unsigned int i=0;i!=Region._StatDataList._InvolvedFluidNum;++i)
		{
          icount++;

          Region._SPHAm[i]._RowID=icount;
          Region._SPHAm[i]._ColID=icount;
          Region._SPHAm[i]._Ele=0.0;

          Region._SPHbu[i]=0.0;
          Region._SPHbv[i]=0.0;

          if(TimeSteps==Region._ControlSPH._StartStep)//计算的初始时间步，将速度置0
			{
              Region._SPHu[i]=0.0;
              Region._SPHv[i]=0.0;
			}
          else//非初始时间步，速度值置为上一时间步的值
			{
              unsigned int ID;
              ID=Region._CalList[i];

              Region._SPHu[i]=Region._PtList[ID]._u;
              Region._SPHv[i]=Region._PtList[ID]._v;
			}
		}
	}
}
