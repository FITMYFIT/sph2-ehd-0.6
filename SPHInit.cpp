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

  //��ز�����0
  Region._ControlSPH._CSPMIflag1=0;//����Ƿ������CSPM���������ķ�ĸ��
  Region._ControlSPH._CSPMIflag2=0;//����Ƿ������CSPM���������ݶȵ�ϵ����


  if(Region._ControlSPH._RunMod==1)//��ʽ�㷨
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

          BasePtPtr->_AccBndx=0.0;//�߽����Ӷ��������Ӳ����ļ��ٶ�
          BasePtPtr->_AccBndy=0.0;//�߽����Ӷ��������Ӳ����ļ��ٶ�

          //BasePtPtr->_NumNegbor=0;//ÿ������֧�����ڵ���������0
          BasePtPtr->_r0=0.0;//��ֹ���룬����particle shift
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

  //��ʽ�㷨
  if(Region._ControlSPH._RunMod==0)
	{
      //�Լ������õ��ı����ó�ֵ
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

          BasePtPtr->_AccBndx=0.0;//�߽����Ӷ��������Ӳ����ļ��ٶ�
          BasePtPtr->_AccBndy=0.0;//�߽����Ӷ��������Ӳ����ļ��ٶ�

          // BasePtPtr->_NumNegbor=0;//ÿ������֧�����ڵ���������0
          BasePtPtr->_r0=0.0;//��ֹ���룬����particle shift
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

      Region._SPHAm.resize(Region._StatDataList._InvolvedFluidNum);//ֻ���������Խ�Ԫ����������Ҫpush_back
      Region._SPHu.resize(Region._StatDataList._InvolvedFluidNum);
      Region._SPHv.resize(Region._StatDataList._InvolvedFluidNum);
      Region._SPHbu.resize(Region._StatDataList._InvolvedFluidNum);
      Region._SPHbv.resize(Region._StatDataList._InvolvedFluidNum);

      //�������е�Ԫ���ú��ʵĳ�ֵ
      unsigned int icount=0;
      for(unsigned int i=0;i!=Region._StatDataList._InvolvedFluidNum;++i)
		{
          icount++;

          Region._SPHAm[i]._RowID=icount;
          Region._SPHAm[i]._ColID=icount;
          Region._SPHAm[i]._Ele=0.0;

          Region._SPHbu[i]=0.0;
          Region._SPHbv[i]=0.0;

          if(TimeSteps==Region._ControlSPH._StartStep)//����ĳ�ʼʱ�䲽�����ٶ���0
			{
              Region._SPHu[i]=0.0;
              Region._SPHv[i]=0.0;
			}
          else//�ǳ�ʼʱ�䲽���ٶ�ֵ��Ϊ��һʱ�䲽��ֵ
			{
              unsigned int ID;
              ID=Region._CalList[i];

              Region._SPHu[i]=Region._PtList[ID]._u;
              Region._SPHv[i]=Region._PtList[ID]._v;
			}
		}
	}
}
