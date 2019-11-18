#include "UpdatePosition.h"

CUpdatePosition::CUpdatePosition()
{
}

CUpdatePosition::~CUpdatePosition()
{
}

void CUpdatePosition::LeapFrogUpdate(CRegion & Region, double DeltaT,double DeltaT1,double TimeSteps)
{
	CBasePt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;
	CBasePt *PtPtr;
	CBasePt * BasePtPtr;
	double DeltaTHf;//DeltT1*0.5
	double DeltaTAver;//(DeltT+DeltT1)*0.5
	unsigned int i,j;
	unsigned int ID;

	double AverRho;
	double uij,vij;
	double XSPHEpsilon;

	//以下变量用于particle shift
	double Maxv;
	double vi;//i粒子的速度模值
	double rij;
	double xij,yij,zij;
	double PtShftBeta;

	DeltaTHf=DeltaT1*0.5;
  DeltaTAver=(DeltaT+DeltaT1)*0.5;

	if(TimeSteps==Region._ControlSPH._StartStep)
    {
      for (i=0;i<Region._PtList.size();i++)
        {
          BasePtPtr=&Region._PtList[i];
          if(BasePtPtr->_Type==enSPHPt)//坐标留待XSPH法计算出速度后再更新
            {
              BasePtPtr->_rhoHf=BasePtPtr->_rho+BasePtPtr->_drho*DeltaTHf;
              BasePtPtr->_uHf=BasePtPtr->_u+BasePtPtr->_du*DeltaTHf;
              BasePtPtr->_vHf=BasePtPtr->_v+BasePtPtr->_dv*DeltaTHf;
              //BasePtPtr->_eHf=BasePtPtr->_e+BasePtPtr->_de*DeltaTHf;
              //BasePtPtr->_THf=BasePtPtr->_T+BasePtPtr->_dT*DeltaTHf;
              BasePtPtr->_hHf=BasePtPtr->_h+BasePtPtr->_dh*DeltaTHf;
            }

          //else if(BasePtPtr->_Type==enNULLPt)//直接更新坐标
          //{
          //	BasePtPtr->_x+=BasePtPtr->_u*DeltaT1;
          //	BasePtPtr->_y+=BasePtPtr->_v*DeltaT1;
          //}

          //if (BasePtPtr->_Type==enDumPt)////dummy粒子，更新密度
          //{
          //	BasePtPtr->_rhoHf=BasePtPtr->_rho+BasePtPtr->_drho*DeltaTHf;
          //}
        }

      if(Region._ControlSPH._XSPHEpsilon!=0)//求XSPH法更新坐标的修正因子，XSPH只用于SPH粒子对之间,刘虎硕士论文，equ.(2.52)
        {

          for(j=0;j<Region._PtList.size();j++)//XSPH-0 XSPH法修正因子置0
            {
              Region._PtList[j]._uXSPHCoef=0;
              Region._PtList[j]._vXSPHCoef=0;
            }

          for(j=0;j<Region._PtPairList.size();j++)//XSPH-1 求各粒子的XSPH修正因子
            {
              if(Region._PtPairList[i]._Type==enSPHPtPair||Region._PtPairList[i]._Type==enSPHNULLPtPair)
                {
                  PtiPtr=Region._PtPairList[j]._PtiPtr;
                  PtjPtr=Region._PtPairList[j]._PtjPtr;

                  KnlPtr=&Region._KnlList[j];

                  uij=PtiPtr->_uHf-PtjPtr->_uHf;
                  vij=PtiPtr->_vHf-PtjPtr->_vHf;

                  AverRho=0.5*(PtiPtr->_rho+PtjPtr->_rho);

                  PtiPtr->_uXSPHCoef+=PtjPtr->_m*uij*KnlPtr->_W/AverRho;
                  PtiPtr->_vXSPHCoef+=PtjPtr->_m*vij*KnlPtr->_W/AverRho;

                  if(PtiPtr!=PtjPtr)
                    {
                      PtjPtr->_uXSPHCoef-=PtiPtr->_m*uij*KnlPtr->_W/AverRho;
                      PtjPtr->_vXSPHCoef-=PtiPtr->_m*vij*KnlPtr->_W/AverRho;
                    }
                }
            }

          XSPHEpsilon=Region._ControlSPH._XSPHEpsilon;

          for(j=0;j<Region._PtList.size();j++)//XSPH-2 XSPH更新坐标
            {
              BasePtPtr=&Region._PtList[j];

              if(BasePtPtr->_Type==enSPHPt)
                {
                  BasePtPtr->_x+=(BasePtPtr->_uHf-XSPHEpsilon*BasePtPtr->_uXSPHCoef)*DeltaT1;
                  BasePtPtr->_y+=(BasePtPtr->_vHf-XSPHEpsilon*BasePtPtr->_vXSPHCoef)*DeltaT1;
                }
            }
        }

      else                                            //不用XSPH法更新enSectionSPH部分的粒子坐标
        {
          for(j=0;j<Region._PtList.size();j++)//非XSPH法更新坐标
            {
              BasePtPtr=&Region._PtList[j];

              if(BasePtPtr->_Type==enSPHPt)
                {
                  BasePtPtr->_x+=BasePtPtr->_uHf*DeltaT1;
                  BasePtPtr->_y+=BasePtPtr->_vHf*DeltaT1;
                }
            }
        }
    }

	else //非第一个时间步
    {
      for (i=0;i<Region._PtList.size();i++)
        {
          BasePtPtr=&Region._PtList[i];

          if(BasePtPtr->_Type==enSPHPt)//对enSectionSPH，只更新速度，不更新坐标，坐标需要判断是否用XSPH法更新
            {
              BasePtPtr->_rhoHf+=BasePtPtr->_drho*DeltaTAver;
              BasePtPtr->_uHf+=BasePtPtr->_du*DeltaTAver;
              BasePtPtr->_vHf+=BasePtPtr->_dv*DeltaTAver;
              //BasePtPtr->_eHf+=BasePtPtr->_de*DeltaTAver;
              //BasePtPtr->_THf+=BasePtPtr->_dT*DeltaTAver;
              BasePtPtr->_hHf+=BasePtPtr->_dh*DeltaTAver;
            }

          //else if(BasePtPtr->_Type==enNULLPt)
          //{
          //	BasePtPtr->_x+=BasePtPtr->_u*DeltaT1;
          //	BasePtPtr->_y+=BasePtPtr->_v*DeltaT1;
          //}

          //if (BasePtPtr->_Type==enDumPt)////dummy粒子，更新密度
          //{
          //	BasePtPtr->_rhoHf=BasePtPtr->_rho+BasePtPtr->_drho*DeltaTHf;
          //}

        }

      if(Region._ControlSPH._XSPHEpsilon!=0)//求XSPH法更新坐标的修正因子，XSPH只用于SPH粒子对及SPH-NULL粒子对之间,刘虎硕士论文，equ.(2.52)
        {

          for(j=0;j<Region._PtList.size();j++)//XSPH-0 XSPH法修正因子置0
            {
              Region._PtList[j]._uXSPHCoef=0;
              Region._PtList[j]._vXSPHCoef=0;
            }

          for(j=0;j<Region._PtPairList.size();j++)//XSPH-1 求各粒子的XSPH修正因子
            {
              if(Region._PtPairList[i]._Type==enSPHPtPair||Region._PtPairList[i]._Type==enSPHNULLPtPair)
                {
                  PtiPtr=Region._PtPairList[j]._PtiPtr;
                  PtjPtr=Region._PtPairList[j]._PtjPtr;

                  KnlPtr=&Region._KnlList[j];

                  uij=PtiPtr->_uHf-PtjPtr->_uHf;
                  vij=PtiPtr->_vHf-PtjPtr->_vHf;

                  AverRho=0.5*(PtiPtr->_rho+PtjPtr->_rho);

                  PtiPtr->_uXSPHCoef+=PtjPtr->_m*uij*KnlPtr->_W/AverRho;
                  PtiPtr->_vXSPHCoef+=PtjPtr->_m*vij*KnlPtr->_W/AverRho;

                  if(PtiPtr!=PtjPtr)
                    {
                      PtjPtr->_uXSPHCoef-=PtiPtr->_m*uij*KnlPtr->_W/AverRho;
                      PtjPtr->_vXSPHCoef-=PtiPtr->_m*vij*KnlPtr->_W/AverRho;
                    }
                }
            }

          XSPHEpsilon=Region._ControlSPH._XSPHEpsilon;

          for(j=0;j<Region._PtList.size();j++)//XSPH-2 XSPH更新坐标
            {
              BasePtPtr=&Region._PtList[j];

              if(BasePtPtr->_Type==enSPHPt)
                {
                  BasePtPtr->_x+=(BasePtPtr->_uHf-XSPHEpsilon*BasePtPtr->_uXSPHCoef)*DeltaT1;
                  BasePtPtr->_y+=(BasePtPtr->_vHf-XSPHEpsilon*BasePtPtr->_vXSPHCoef)*DeltaT1;
                }
            }
        }

      //particle shift algrithm,2012_A robust weakly compressible SPH method and its comparison with an incompressible SPH
      //deltai=beta*SUM(xij/rij^3*r0^2*vmax)
      else if (Region._ControlSPH._PtShftBeta!=0)
        {
          //1.找到最大速度
          Maxv=0.0;
          for (i=0;i!=Region._CalList.size();++i)
            {
              ID=Region._CalList[i];

              vi=sqrt(Region._PtList[ID]._uHf*Region._PtList[ID]._uHf+Region._PtList[ID]._vHf*Region._PtList[ID]._vHf);

              if (vi>Maxv)
                {
                  Maxv=vi;
                }
            }

          //2.计算截止距离（cutoff distance）
          //r0=sum(rij/N)
          for (i=0;i!=Region._PtPairList.size();++i)
            {
              PtiPtr=Region._PtPairList[i]._PtiPtr;
              PtjPtr=Region._PtPairList[i]._PtjPtr;

              rij=sqrt(Region._PtPairList[i]._driac2);

              PtiPtr->_r0+=rij/PtiPtr->_NumNegbor;

              if (PtiPtr!=PtjPtr)
                {
                  PtjPtr->_r0+=rij/PtjPtr->_NumNegbor;
                }
            }

          //3.计算修正系数
          PtShftBeta=Region._ControlSPH._PtShftBeta;
          for (i=0;i!=Region._PtPairList.size();++i)
            {
              PtiPtr=Region._PtPairList[i]._PtiPtr;
              PtjPtr=Region._PtPairList[i]._PtjPtr;

              if (PtiPtr!=PtjPtr)
                {
                  rij=sqrt(Region._PtPairList[i]._driac2);
                  xij=Region._PtPairList[i]._xij;          
                  yij=Region._PtPairList[i]._yij;

                  PtiPtr->_PtShftCoefu+=PtShftBeta*xij/(rij*rij*rij)*PtiPtr->_r0*PtiPtr->_r0*Maxv;
                  PtiPtr->_PtShftCoefv+=PtShftBeta*yij/(rij*rij*rij)*PtiPtr->_r0*PtiPtr->_r0*Maxv;

                  PtjPtr->_PtShftCoefu-=PtShftBeta*xij/(rij*rij*rij)*PtjPtr->_r0*PtjPtr->_r0*Maxv;
                  PtjPtr->_PtShftCoefv-=PtShftBeta*yij/(rij*rij*rij)*PtjPtr->_r0*PtjPtr->_r0*Maxv;
                }
            }

          for(j=0;j<Region._PtList.size();j++)//Particle Shift更新坐标
            {
              BasePtPtr=&Region._PtList[j];

              if(BasePtPtr->_Type==enSPHPt)
                {
                  BasePtPtr->_x+=(BasePtPtr->_uHf+BasePtPtr->_PtShftCoefu)*DeltaT1;
                  BasePtPtr->_y+=(BasePtPtr->_vHf+BasePtPtr->_PtShftCoefv)*DeltaT1;
                }
            }
        }
      //particle shift algrithm,2012_A robust weakly compressible SPH method and its comparison with an incompressible SPH
      //deltai=beta*SUM(xij/rij^3*r0^2*vmax)
      else if (Region._ControlSPH._PtShftBeta!=0)
        {
          //1.找到最大速度
          Maxv=0.0;
          for (i=0;i!=Region._CalList.size();++i)
            {
              ID=Region._CalList[i];

              vi=sqrt(Region._PtList[ID]._uHf*Region._PtList[ID]._uHf+Region._PtList[ID]._vHf*Region._PtList[ID]._vHf);

              if (vi>Maxv)
                {
                  Maxv=vi;
                }
            }

          //2.计算截止距离（cutoff distance）
          //r0=sum(rij/N)
          for (i=0;i!=Region._PtPairList.size();++i)
            {
              PtiPtr=Region._PtPairList[i]._PtiPtr;
              PtjPtr=Region._PtPairList[i]._PtjPtr;

              rij=sqrt(Region._PtPairList[i]._driac2);

              PtiPtr->_r0+=rij/PtiPtr->_NumNegbor;

              if (PtiPtr!=PtjPtr)
                {
                  PtjPtr->_r0+=rij/PtjPtr->_NumNegbor;
                }
            }

          //3.计算修正系数
          PtShftBeta=Region._ControlSPH._PtShftBeta;
          for (i=0;i!=Region._PtPairList.size();++i)
            {
              PtiPtr=Region._PtPairList[i]._PtiPtr;
              PtjPtr=Region._PtPairList[i]._PtjPtr;

              if (PtiPtr!=PtjPtr)
                {
                  rij=sqrt(Region._PtPairList[i]._driac2);
                  xij=Region._PtPairList[i]._xij;
                  yij=Region._PtPairList[i]._yij;

                  PtiPtr->_PtShftCoefu+=PtShftBeta*xij/(rij*rij*rij)*PtiPtr->_r0*PtiPtr->_r0*Maxv;
                  PtiPtr->_PtShftCoefv+=PtShftBeta*yij/(rij*rij*rij)*PtiPtr->_r0*PtiPtr->_r0*Maxv;

                  PtjPtr->_PtShftCoefu-=PtShftBeta*xij/(rij*rij*rij)*PtjPtr->_r0*PtjPtr->_r0*Maxv;
                  PtjPtr->_PtShftCoefv-=PtShftBeta*yij/(rij*rij*rij)*PtjPtr->_r0*PtjPtr->_r0*Maxv;
                }
            }

          for(j=0;j<Region._PtList.size();j++)//Particle Shift更新坐标
            {
              BasePtPtr=&Region._PtList[j];

              if(BasePtPtr->_Type==enSPHPt)
                {
                  BasePtPtr->_x+=(BasePtPtr->_uHf+BasePtPtr->_PtShftCoefu)*DeltaT1;
                  BasePtPtr->_y+=(BasePtPtr->_vHf+BasePtPtr->_PtShftCoefv)*DeltaT1;
                }
            }
        }

      else                                            //不用XSPH法更新enSectionSPH部分的粒子坐标
        {
          for(j=0;j<Region._PtList.size();j++)//非XSPH法更新坐标
            {
              BasePtPtr=&Region._PtList[j];

              if(BasePtPtr->_Type==enSPHPt)
                {
                  BasePtPtr->_x+=BasePtPtr->_uHf*DeltaT1;
                  BasePtPtr->_y+=BasePtPtr->_vHf*DeltaT1;
                }
            }
        }
    }

	for (i=0;i<Region._PtList.size();i++)
    {
      BasePtPtr=&Region._PtList[i];
      if(BasePtPtr->_Type==enSPHPt)
        {
          BasePtPtr->_rho=BasePtPtr->_rhoHf+BasePtPtr->_drho*DeltaTHf;
          BasePtPtr->_u=BasePtPtr->_uHf+BasePtPtr->_du*DeltaTHf;
          BasePtPtr->_v=BasePtPtr->_vHf+BasePtPtr->_dv*DeltaTHf;
          //BasePtPtr->_e=BasePtPtr->_eHf+BasePtPtr->_de*DeltaTHf;
          //BasePtPtr->_T=BasePtPtr->_THf+BasePtPtr->_dT*DeltaTHf;
          BasePtPtr->_h=BasePtPtr->_hHf+BasePtPtr->_dh*DeltaTHf;
          BasePtPtr->_r=2.0*BasePtPtr->_h;
        }

      //if (BasePtPtr->_Type==enNULLPt)//更新Null粒子的密度，以获得压力
      //{
      //	if (BasePtPtr->_PID!=3)//运动的上顶没有压力
      //	{
      //		BasePtPtr->_rho+=BasePtPtr->_drho*Region._ControlSPH._DeltaT;
      //	}
      //}

      //if (BasePtPtr->_Type==enDumPt)////dummy粒子，更新密度
      //{
      //	BasePtPtr->_rho=BasePtPtr->_rhoHf+BasePtPtr->_drho*DeltaTHf;
      //}

    }

	//周期性边界，如果某个粒子超出了周期性边界，则移动该粒子
	if (Region._ControlSPH._PerdBnd==1)
    {
      double MAXX,MINX;//用于周期性边界条件
      double dp;//粒子间距

      dp=sqrt(Region._PtList[0]._Volume);

      MAXX=Region._ControlSPH._PerdBndMaxX;
      MINX=Region._ControlSPH._PerdBndMinX;

      for (i=0;i!=Region._PtList.size();++i)
        {
          if (Region._PtList[i]._x>=MAXX+dp)
            {
              Region._PtList[i]._x-=(MAXX-MINX);
            }
        }
    }



	////dam break ETSIN 2014,挡板运动
	//double Vgate;//挡板运动速度
	//for (unsigned int i=0;i!=Region._PtList.size();++i)
	//{
	//	BasePtPtr=&Region._PtList[i];

	//	if(BasePtPtr->_PID==2)//挡板的PID=2
	//	{
	//		Vgate=BasePtPtr->_v;

	//		break;
	//	}
	//}

	////if(TimeSteps*DeltaT<Region._ControlSPH._TDamp)//考虑damp过程
	//{
	//	if (TimeSteps*DeltaT*Vgate<0.63)
	//	{
	//		for (i=0;i!=Region._PtList.size();++i)
	//		{
	//			BasePtPtr=&Region._PtList[i];

	//			if(BasePtPtr->_PID==2)//挡板的PID=2
	//			{
	//				BasePtPtr->_y+=BasePtPtr->_v*Region._ControlSPH._DeltaT;
	//			}
	//		}
	//	}
	//}

	////control disk算例，part 2做匀速直线运动
	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	PtPtr=&Region._PtList[i];
	//	if (PtPtr->_PID==2)
	//	{
	//		PtPtr->_x+=PtPtr->_u*Region._ControlSPH._DeltaT;
	//		PtPtr->_y+=PtPtr->_v*Region._ControlSPH._DeltaT;
	//	}

	//}




	////盒子在流场中运动算例，盒子以a=1m/s2的加速度运动至1m/s，而后做匀速运动
	////盒子是part2
	////1.找出现在盒子的运动速度
	//double unow;
	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	if (Region._PtList[i]._PID==2)
	//	{
	//		unow=Region._PtList[i]._u;
	//		break;
	//	}
	//}

	////2.根据运动速度判断盒子处于加速或是匀速运动状态
	//double acc;//加速度
	//double deltat;//时间步长
	//deltat=Region._ControlSPH._DeltaT;
	//acc=1.0;

	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	PtPtr=&Region._PtList[i];
	//	if (PtPtr->_PID==2)
	//	{
	//		if (PtPtr->_u<PtPtr->_vwall)//处于加速运动状态
	//		{
	//			PtPtr->_x+=(PtPtr->_u+0.5*acc*deltat)*deltat;
	//			PtPtr->_u+=acc*deltat;
	//		}

	//		else//处于匀速运动状态
	//		{
	//			PtPtr->_x+=PtPtr->_u*deltat;
	//		}
	//	}
	//}
}

void CUpdatePosition::ImplicitUpdate(CRegion & Region)
{
	int i,j;
	unsigned int ID,ID2;
	unsigned int IDi,IDj;
	CBasePt *PtPtr;
	CBasePt * PtiPtr,*PtjPtr;
	CKnl * KnlPtr;
	double uij,vij;
	double xx;
	double accuracy;
	double AverRho;
	float XSPHEpsilon;

	//非XSPH法推进位置坐标
	if(Region._ControlSPH._XSPHEpsilon==0)
	{
		for(i=0;i<Region._PtList.size();i++)
		{
			PtPtr=&Region._PtList[i];

			if(PtPtr->_Type==enSPHPt&&PtPtr->_Iflag==1)
			{
				ID=PtPtr->_ID2;

				PtPtr->_x+=0.5*Region._ControlSPH._DeltaT*(PtPtr->_u+Region._SPHu[ID-1]);
				PtPtr->_y+=0.5*Region._ControlSPH._DeltaT*(PtPtr->_v+Region._SPHv[ID-1]);
			}
			else//Null粒子
			{
				//PtPtr->_x+=Region._ControlSPH._DeltaT*PtPtr->_u;
				//PtPtr->_y+=Region._ControlSPH._DeltaT*PtPtr->_v;
			}
		}
	}

	//XSPH法推进粒子坐标
	else
	{
		for(j=0;j<Region._PtList.size();j++)//XSPH-0 XSPH法修正因子置0
		{
			if((Region._PtList[j]._Type==enSPHPt||Region._PtList[j]._Type==enNULLPt)&&Region._PtList[j]._Iflag==1)
			{
				Region._PtList[j]._uXSPHCoef=0.0;
				Region._PtList[j]._vXSPHCoef=0.0;
			}
		}

		for(j=0;j<Region._PtPairList.size();j++)//XSPH-1 求各粒子的XSPH修正因子
		{
			if(Region._PtPairList[j]._Type!=enSPHBndPtPair)
			{
				PtiPtr=Region._PtPairList[j]._PtiPtr;
				PtjPtr=Region._PtPairList[j]._PtjPtr;

				KnlPtr=&Region._KnlList[j];

				IDi=PtiPtr->_ID2;
				IDj=PtjPtr->_ID2;

				uij=0.5*(PtiPtr->_u+Region._SPHu[IDi-1])-0.5*(PtjPtr->_u+Region._SPHu[IDj-1]);//用前后时间步的速度的平均值
				vij=0.5*(PtiPtr->_v+Region._SPHv[IDi-1])-0.5*(PtjPtr->_v+Region._SPHv[IDj-1]);

				AverRho=0.5*(PtiPtr->_rho+PtjPtr->_rho);

				PtiPtr->_uXSPHCoef+=PtjPtr->_m*uij*KnlPtr->_W/AverRho;
				PtiPtr->_vXSPHCoef+=PtjPtr->_m*vij*KnlPtr->_W/AverRho;

				if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)
				{
					PtjPtr->_uXSPHCoef-=PtiPtr->_m*uij*KnlPtr->_W/AverRho;
					PtjPtr->_vXSPHCoef-=PtiPtr->_m*vij*KnlPtr->_W/AverRho;
				}
			}
		}

		XSPHEpsilon=Region._ControlSPH._XSPHEpsilon;

		for(j=0;j<Region._PtList.size();j++)//XSPH-2 XSPH更新坐标
		{
			PtPtr=&Region._PtList[j];

			if(Region._PtList[j]._Type==enSPHPt&&Region._PtList[j]._Iflag==1)
			{
				ID=PtPtr->_ID2;

				PtPtr->_x+=(0.5*(PtPtr->_u+Region._SPHu[ID-1])-XSPHEpsilon*PtPtr->_uXSPHCoef)*Region._ControlSPH._DeltaT;
				PtPtr->_y+=(0.5*(PtPtr->_v+Region._SPHv[ID-1])-XSPHEpsilon*PtPtr->_vXSPHCoef)*Region._ControlSPH._DeltaT;
			}

			else//Null粒子及Bnd粒子
			{
				//PtPtr->_x+=Region._ControlSPH._DeltaT*PtPtr->_u;
				//PtPtr->_y+=Region._ControlSPH._DeltaT*PtPtr->_v;
			}
		}
	}

	//将计算出的速度值转存回SPHPtList中
	for(i=0;i<Region._CalList.size();i++)
	{
		ID=Region._CalList[i];

		Region._PtList[ID]._u=Region._SPHu[i];
		Region._PtList[ID]._v=Region._SPHv[i];
	}

	//周期性边界，如果某个粒子超出了周期性边界，则移动该粒子
	if (Region._ControlSPH._PerdBnd==1)
	{
		double MAXX,MINX;//用于周期性边界条件
		double dp;//粒子间距

		dp=sqrt(Region._PtList[0]._Volume);

		MAXX=Region._ControlSPH._PerdBndMaxX;
		MINX=Region._ControlSPH._PerdBndMinX;

		for (i=0;i!=Region._PtList.size();++i)
		{
			if (Region._PtList[i]._x>=MAXX+dp)
			{
				Region._PtList[i]._x-=(MAXX-MINX);
			}
		}
	}


	////盒子在流场中运动算例，盒子以a=1m/s2的加速度运动至1m/s，而后做匀速运动
	////盒子是part2
	////1.找出现在盒子的运动速度
	//double unow;
	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	if (Region._PtList[i]._PID==2)
	//	{
	//		unow=Region._PtList[i]._u;
	//		break;
	//	}
	//}

	////2.根据运动速度判断盒子处于加速或是匀速运动状态
	//double acc;//加速度
	//double deltat;//时间步长
	//deltat=Region._ControlSPH._DeltaT;
	//acc=1.0;

	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	PtPtr=&Region._PtList[i];
	//	if (PtPtr->_PID==2)
	//	{
	//		if (PtPtr->_u<1.0)//处于加速运动状态
	//		{
	//			PtPtr->_x+=(PtPtr->_u+0.5*acc*deltat)*deltat;
	//			PtPtr->_u+=acc*deltat;
	//		}

	//		else//处于匀速运动状态
	//		{
	//			PtPtr->_x+=PtPtr->_u*deltat;
	//		}
	//	}
	//}
}

