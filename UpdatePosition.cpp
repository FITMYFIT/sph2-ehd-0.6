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

	//���±�������particle shift
	double Maxv;
	double vi;//i���ӵ��ٶ�ģֵ
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
          if(BasePtPtr->_Type==enSPHPt)//��������XSPH��������ٶȺ��ٸ���
            {
              BasePtPtr->_rhoHf=BasePtPtr->_rho+BasePtPtr->_drho*DeltaTHf;
              BasePtPtr->_uHf=BasePtPtr->_u+BasePtPtr->_du*DeltaTHf;
              BasePtPtr->_vHf=BasePtPtr->_v+BasePtPtr->_dv*DeltaTHf;
              //BasePtPtr->_eHf=BasePtPtr->_e+BasePtPtr->_de*DeltaTHf;
              //BasePtPtr->_THf=BasePtPtr->_T+BasePtPtr->_dT*DeltaTHf;
              BasePtPtr->_hHf=BasePtPtr->_h+BasePtPtr->_dh*DeltaTHf;
            }

          //else if(BasePtPtr->_Type==enNULLPt)//ֱ�Ӹ�������
          //{
          //	BasePtPtr->_x+=BasePtPtr->_u*DeltaT1;
          //	BasePtPtr->_y+=BasePtPtr->_v*DeltaT1;
          //}

          //if (BasePtPtr->_Type==enDumPt)////dummy���ӣ������ܶ�
          //{
          //	BasePtPtr->_rhoHf=BasePtPtr->_rho+BasePtPtr->_drho*DeltaTHf;
          //}
        }

      if(Region._ControlSPH._XSPHEpsilon!=0)//��XSPH������������������ӣ�XSPHֻ����SPH���Ӷ�֮��,����˶ʿ���ģ�equ.(2.52)
        {

          for(j=0;j<Region._PtList.size();j++)//XSPH-0 XSPH������������0
            {
              Region._PtList[j]._uXSPHCoef=0;
              Region._PtList[j]._vXSPHCoef=0;
            }

          for(j=0;j<Region._PtPairList.size();j++)//XSPH-1 ������ӵ�XSPH��������
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

          for(j=0;j<Region._PtList.size();j++)//XSPH-2 XSPH��������
            {
              BasePtPtr=&Region._PtList[j];

              if(BasePtPtr->_Type==enSPHPt)
                {
                  BasePtPtr->_x+=(BasePtPtr->_uHf-XSPHEpsilon*BasePtPtr->_uXSPHCoef)*DeltaT1;
                  BasePtPtr->_y+=(BasePtPtr->_vHf-XSPHEpsilon*BasePtPtr->_vXSPHCoef)*DeltaT1;
                }
            }
        }

      else                                            //����XSPH������enSectionSPH���ֵ���������
        {
          for(j=0;j<Region._PtList.size();j++)//��XSPH����������
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

	else //�ǵ�һ��ʱ�䲽
    {
      for (i=0;i<Region._PtList.size();i++)
        {
          BasePtPtr=&Region._PtList[i];

          if(BasePtPtr->_Type==enSPHPt)//��enSectionSPH��ֻ�����ٶȣ����������꣬������Ҫ�ж��Ƿ���XSPH������
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

          //if (BasePtPtr->_Type==enDumPt)////dummy���ӣ������ܶ�
          //{
          //	BasePtPtr->_rhoHf=BasePtPtr->_rho+BasePtPtr->_drho*DeltaTHf;
          //}

        }

      if(Region._ControlSPH._XSPHEpsilon!=0)//��XSPH������������������ӣ�XSPHֻ����SPH���ӶԼ�SPH-NULL���Ӷ�֮��,����˶ʿ���ģ�equ.(2.52)
        {

          for(j=0;j<Region._PtList.size();j++)//XSPH-0 XSPH������������0
            {
              Region._PtList[j]._uXSPHCoef=0;
              Region._PtList[j]._vXSPHCoef=0;
            }

          for(j=0;j<Region._PtPairList.size();j++)//XSPH-1 ������ӵ�XSPH��������
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

          for(j=0;j<Region._PtList.size();j++)//XSPH-2 XSPH��������
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
          //1.�ҵ�����ٶ�
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

          //2.�����ֹ���루cutoff distance��
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

          //3.��������ϵ��
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

          for(j=0;j<Region._PtList.size();j++)//Particle Shift��������
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
          //1.�ҵ�����ٶ�
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

          //2.�����ֹ���루cutoff distance��
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

          //3.��������ϵ��
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

          for(j=0;j<Region._PtList.size();j++)//Particle Shift��������
            {
              BasePtPtr=&Region._PtList[j];

              if(BasePtPtr->_Type==enSPHPt)
                {
                  BasePtPtr->_x+=(BasePtPtr->_uHf+BasePtPtr->_PtShftCoefu)*DeltaT1;
                  BasePtPtr->_y+=(BasePtPtr->_vHf+BasePtPtr->_PtShftCoefv)*DeltaT1;
                }
            }
        }

      else                                            //����XSPH������enSectionSPH���ֵ���������
        {
          for(j=0;j<Region._PtList.size();j++)//��XSPH����������
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

      //if (BasePtPtr->_Type==enNULLPt)//����Null���ӵ��ܶȣ��Ի��ѹ��
      //{
      //	if (BasePtPtr->_PID!=3)//�˶����϶�û��ѹ��
      //	{
      //		BasePtPtr->_rho+=BasePtPtr->_drho*Region._ControlSPH._DeltaT;
      //	}
      //}

      //if (BasePtPtr->_Type==enDumPt)////dummy���ӣ������ܶ�
      //{
      //	BasePtPtr->_rho=BasePtPtr->_rhoHf+BasePtPtr->_drho*DeltaTHf;
      //}

    }

	//�����Ա߽磬���ĳ�����ӳ����������Ա߽磬���ƶ�������
	if (Region._ControlSPH._PerdBnd==1)
    {
      double MAXX,MINX;//���������Ա߽�����
      double dp;//���Ӽ��

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



	////dam break ETSIN 2014,�����˶�
	//double Vgate;//�����˶��ٶ�
	//for (unsigned int i=0;i!=Region._PtList.size();++i)
	//{
	//	BasePtPtr=&Region._PtList[i];

	//	if(BasePtPtr->_PID==2)//�����PID=2
	//	{
	//		Vgate=BasePtPtr->_v;

	//		break;
	//	}
	//}

	////if(TimeSteps*DeltaT<Region._ControlSPH._TDamp)//����damp����
	//{
	//	if (TimeSteps*DeltaT*Vgate<0.63)
	//	{
	//		for (i=0;i!=Region._PtList.size();++i)
	//		{
	//			BasePtPtr=&Region._PtList[i];

	//			if(BasePtPtr->_PID==2)//�����PID=2
	//			{
	//				BasePtPtr->_y+=BasePtPtr->_v*Region._ControlSPH._DeltaT;
	//			}
	//		}
	//	}
	//}

	////control disk������part 2������ֱ���˶�
	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	PtPtr=&Region._PtList[i];
	//	if (PtPtr->_PID==2)
	//	{
	//		PtPtr->_x+=PtPtr->_u*Region._ControlSPH._DeltaT;
	//		PtPtr->_y+=PtPtr->_v*Region._ControlSPH._DeltaT;
	//	}

	//}




	////�������������˶�������������a=1m/s2�ļ��ٶ��˶���1m/s�������������˶�
	////������part2
	////1.�ҳ����ں��ӵ��˶��ٶ�
	//double unow;
	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	if (Region._PtList[i]._PID==2)
	//	{
	//		unow=Region._PtList[i]._u;
	//		break;
	//	}
	//}

	////2.�����˶��ٶ��жϺ��Ӵ��ڼ��ٻ��������˶�״̬
	//double acc;//���ٶ�
	//double deltat;//ʱ�䲽��
	//deltat=Region._ControlSPH._DeltaT;
	//acc=1.0;

	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	PtPtr=&Region._PtList[i];
	//	if (PtPtr->_PID==2)
	//	{
	//		if (PtPtr->_u<PtPtr->_vwall)//���ڼ����˶�״̬
	//		{
	//			PtPtr->_x+=(PtPtr->_u+0.5*acc*deltat)*deltat;
	//			PtPtr->_u+=acc*deltat;
	//		}

	//		else//���������˶�״̬
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

	//��XSPH���ƽ�λ������
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
			else//Null����
			{
				//PtPtr->_x+=Region._ControlSPH._DeltaT*PtPtr->_u;
				//PtPtr->_y+=Region._ControlSPH._DeltaT*PtPtr->_v;
			}
		}
	}

	//XSPH���ƽ���������
	else
	{
		for(j=0;j<Region._PtList.size();j++)//XSPH-0 XSPH������������0
		{
			if((Region._PtList[j]._Type==enSPHPt||Region._PtList[j]._Type==enNULLPt)&&Region._PtList[j]._Iflag==1)
			{
				Region._PtList[j]._uXSPHCoef=0.0;
				Region._PtList[j]._vXSPHCoef=0.0;
			}
		}

		for(j=0;j<Region._PtPairList.size();j++)//XSPH-1 ������ӵ�XSPH��������
		{
			if(Region._PtPairList[j]._Type!=enSPHBndPtPair)
			{
				PtiPtr=Region._PtPairList[j]._PtiPtr;
				PtjPtr=Region._PtPairList[j]._PtjPtr;

				KnlPtr=&Region._KnlList[j];

				IDi=PtiPtr->_ID2;
				IDj=PtjPtr->_ID2;

				uij=0.5*(PtiPtr->_u+Region._SPHu[IDi-1])-0.5*(PtjPtr->_u+Region._SPHu[IDj-1]);//��ǰ��ʱ�䲽���ٶȵ�ƽ��ֵ
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

		for(j=0;j<Region._PtList.size();j++)//XSPH-2 XSPH��������
		{
			PtPtr=&Region._PtList[j];

			if(Region._PtList[j]._Type==enSPHPt&&Region._PtList[j]._Iflag==1)
			{
				ID=PtPtr->_ID2;

				PtPtr->_x+=(0.5*(PtPtr->_u+Region._SPHu[ID-1])-XSPHEpsilon*PtPtr->_uXSPHCoef)*Region._ControlSPH._DeltaT;
				PtPtr->_y+=(0.5*(PtPtr->_v+Region._SPHv[ID-1])-XSPHEpsilon*PtPtr->_vXSPHCoef)*Region._ControlSPH._DeltaT;
			}

			else//Null���Ӽ�Bnd����
			{
				//PtPtr->_x+=Region._ControlSPH._DeltaT*PtPtr->_u;
				//PtPtr->_y+=Region._ControlSPH._DeltaT*PtPtr->_v;
			}
		}
	}

	//����������ٶ�ֵת���SPHPtList��
	for(i=0;i<Region._CalList.size();i++)
	{
		ID=Region._CalList[i];

		Region._PtList[ID]._u=Region._SPHu[i];
		Region._PtList[ID]._v=Region._SPHv[i];
	}

	//�����Ա߽磬���ĳ�����ӳ����������Ա߽磬���ƶ�������
	if (Region._ControlSPH._PerdBnd==1)
	{
		double MAXX,MINX;//���������Ա߽�����
		double dp;//���Ӽ��

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


	////�������������˶�������������a=1m/s2�ļ��ٶ��˶���1m/s�������������˶�
	////������part2
	////1.�ҳ����ں��ӵ��˶��ٶ�
	//double unow;
	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	if (Region._PtList[i]._PID==2)
	//	{
	//		unow=Region._PtList[i]._u;
	//		break;
	//	}
	//}

	////2.�����˶��ٶ��жϺ��Ӵ��ڼ��ٻ��������˶�״̬
	//double acc;//���ٶ�
	//double deltat;//ʱ�䲽��
	//deltat=Region._ControlSPH._DeltaT;
	//acc=1.0;

	//for (i=0;i!=Region._PtList.size();++i)
	//{
	//	PtPtr=&Region._PtList[i];
	//	if (PtPtr->_PID==2)
	//	{
	//		if (PtPtr->_u<1.0)//���ڼ����˶�״̬
	//		{
	//			PtPtr->_x+=(PtPtr->_u+0.5*acc*deltat)*deltat;
	//			PtPtr->_u+=acc*deltat;
	//		}

	//		else//���������˶�״̬
	//		{
	//			PtPtr->_x+=PtPtr->_u*deltat;
	//		}
	//	}
	//}
}

