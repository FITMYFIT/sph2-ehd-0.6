#include "ContinuityEqu.h"

CContinuityEqu::CContinuityEqu()
{
}

CContinuityEqu::~CContinuityEqu()
{
}

void CContinuityEqu::Solve(CRegion &Region,unsigned int TimeSteps)
{
	CBasePt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;

	double uij,vij;
	double temp;
	unsigned int IDi,IDj;
	double mi,mj;
	double xij,yij,zij;
	double rhoi,rhoj;
	double ADdelta;
	double wcx,wcy;
	double *gradrhox,*gradrhoy,*gradrhoz;//密度梯度，用于在连续性方程中施加人工耗散项
	double psix,psiy,psiz;//人工耗散项中的psiij的分量
	double hij,cs;//人工耗散项中的光滑长度与声速
	double modulev;//速度的模

	for (size_t i=0;i<Region._PtPairList.size();i++)
	{
		if(Region._PtPairList[i]._Type==enSPHPtPair||Region._PtPairList[i]._Type==enSPHNULLPtPair)
		{
			PtiPtr=Region._PtPairList[i]._PtiPtr;
			PtjPtr=Region._PtPairList[i]._PtjPtr;
			KnlPtr=&Region._KnlList[i];

			uij=PtiPtr->_u-PtjPtr->_u;
			vij=PtiPtr->_v-PtjPtr->_v;

			//temp=uij*KnlPtr->_Wxj+vij*KnlPtr->_Wyj;//用粒子j的光滑长度计算出的核函数信息
			temp=uij*KnlPtr->_Wx+vij*KnlPtr->_Wy;//用粒子j的光滑长度计算出的核函数信息

			PtiPtr->_drho+=PtiPtr->_m*temp;

			if(PtiPtr!=PtjPtr)
			{
				//temp=uij*KnlPtr->_Wx+vij*KnlPtr->_Wy;//用粒子i的光滑长度计算出的核函数信息
				PtjPtr->_drho+=PtjPtr->_m*temp;
			}
		}

		if(Region._PtPairList[i]._Type==enSPHDumPtPair)
		{
			PtiPtr=Region._PtPairList[i]._PtiPtr;
			PtjPtr=Region._PtPairList[i]._PtjPtr;
			KnlPtr=&Region._KnlList[i];

			uij=PtiPtr->_u-PtjPtr->_u;
			vij=PtiPtr->_v-PtjPtr->_v;
			

			modulev=sqrt(PtiPtr->_u*PtiPtr->_u+PtiPtr->_v*PtiPtr->_v);

			//if (fabs(uij*PtjPtr->_Nnx+vij*PtjPtr->_Nny)>=0.00001*modulev)
			{
				//temp=uij*KnlPtr->_Wxj+vij*KnlPtr->_Wyj;//用粒子j的光滑长度计算出的核函数信息
				temp=uij*KnlPtr->_Wx+vij*KnlPtr->_Wy;//用粒子j的光滑长度计算出的核函数信息

				PtiPtr->_drho+=PtiPtr->_m*temp;
			}
		}
	}

	//Artificial Diffusive项（人工耗散项）2011—S Marrone
	if(Region._ControlSPH._ADDelta!=0)
	{
		//1.求密度的梯度（利用修正核函数梯度）
		//Marrone S Equ.(6)
		gradrhox=(double *)malloc(Region._PtList.size()*sizeof(double));
		gradrhoy=(double *)malloc(Region._PtList.size()*sizeof(double));

		for(unsigned int i=0;i<Region._PtList.size();i++)
		{
			gradrhox[i]=0.0;
			gradrhoy[i]=0.0;
		}

		for(unsigned int i=0;i!=Region._PtPairList.size();i++)
		{
			if(Region._PtPairList[i]._Type==enSPHPtPair)
			{
				PtiPtr=Region._PtPairList[i]._PtiPtr;
				PtjPtr=Region._PtPairList[i]._PtjPtr;
				KnlPtr=&Region._KnlList[i];

				IDi=PtiPtr->_ID;
				IDj=PtjPtr->_ID;

				mi=PtiPtr->_m;
				mj=PtjPtr->_m;

				rhoi=PtiPtr->_rho;
				rhoj=PtjPtr->_rho;

				_GetKnlListC.GetCrctKnlGrad(Region,PtiPtr,KnlPtr,wcx,wcy);

				gradrhox[IDi-1]+=mj/rhoj*(rhoj-rhoi)*wcx;
				gradrhoy[IDi-1]+=mj/rhoj*(rhoj-rhoi)*wcy;

				if(PtiPtr!=PtjPtr)
				{
					_GetKnlListC.GetCrctKnlGrad(Region,PtjPtr,KnlPtr,wcx,wcy);

					gradrhox[IDj-1]-=mi/rhoi*(rhoi-rhoj)*wcx;
					gradrhoy[IDj-1]-=mi/rhoi*(rhoi-rhoj)*wcy;
				}
			}
		}

		//2.求人工耗散项
		cs=Region._ControlSPH._Cs;
		ADdelta=Region._ControlSPH._ADDelta;
		for(unsigned int i=0;i!=Region._PtPairList.size();i++)
		{
			if(Region._PtPairList[i]._Type==enSPHPtPair)
			{
				PtiPtr=Region._PtPairList[i]._PtiPtr;
				PtjPtr=Region._PtPairList[i]._PtjPtr;
				KnlPtr=&Region._KnlList[i];

				IDi=PtiPtr->_ID;
				IDj=PtjPtr->_ID;

				mi=PtiPtr->_m;
				mj=PtjPtr->_m;

				rhoi=PtiPtr->_rho;
				rhoj=PtjPtr->_rho;

        xij=Region._PtPairList[i]._xij;
        yij=Region._PtPairList[i]._yij;

	          if(PtiPtr!=PtjPtr)
              {
                psix=-2*(rhoj-rhoi)*xij/(xij*xij+yij*yij)-(gradrhox[IDi-1]+gradrhox[IDj-1]);
                psiy=-2*(rhoj-rhoi)*yij/(xij*xij+yij*yij)-(gradrhoy[IDi-1]+gradrhoy[IDj-1]);

                hij=0.5*(PtiPtr->_h+PtjPtr->_h);

                PtiPtr->_drho+=ADdelta*hij*cs*mj/rhoj*(psix*KnlPtr->_Wx+psiy*KnlPtr->_Wy);

                PtjPtr->_drho-=ADdelta*hij*cs*mi/rhoi*(psix*KnlPtr->_Wx+psiy*KnlPtr->_Wy);
              }
			}
		}

		free(gradrhox);
		free(gradrhoy);
	}

	if(Region._ControlSPH._RunMod==0)//隐式算法，更新粒子的密度，显式算法的密度在LeapFrog处更新
	{
		double DeltaT;
		DeltaT=Region._ControlSPH._DeltaT;
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			if(Region._PtList[i]._Type==enSPHPt&&Region._PtList[i]._Iflag==1)
			{
				Region._PtList[i]._rho+=Region._PtList[i]._drho*DeltaT;
			}

			if(Region._PtList[i]._Type==enNULLPt&&Region._PtList[i]._Iflag==1)//Null粒子的密度也更新，以期获得Null粒子的压力
			{
				Region._PtList[i]._rho+=Region._PtList[i]._drho*DeltaT;
			}
		}
	}





	//按指定的步数重构密度
	if(Region._ControlSPH._DensRenormSteps!=0)
	{
		if(TimeSteps!=Region._ControlSPH._StartStep&&TimeSteps%Region._ControlSPH._DensRenormSteps==0)
		{
			if(Region._ControlSPH._CSPMIflag1==0)//没有计算CSPM修正系数
			{
				_CSPM.GetCSPMFunCorctCoef(Region);
			}

			//CSPM-3 将粒子的密度置0
			for(unsigned int i=0;i<Region._PtList.size();i++)
			{
				if(Region._PtList[i]._Type==enSPHPt&&Region._PtList[i]._Iflag==1)
				{
					Region._PtList[i]._rho=0.0;
				}
			}

			//CSPM-4 插值求密度
			for (unsigned int i=0;i<Region._PtPairList.size();i++)
			{
				if(Region._PtPairList[i]._Type==enSPHPtPair)
				{
					PtiPtr=Region._PtPairList[i]._PtiPtr;
					PtjPtr=Region._PtPairList[i]._PtjPtr;
					KnlPtr=&Region._KnlList[i];

					PtiPtr->_rho+=PtjPtr->_m*KnlPtr->_W;

					if(PtiPtr!=PtjPtr)
					{
						PtjPtr->_rho+=PtiPtr->_m*KnlPtr->_W;
					}
				}

				if(Region._PtPairList[i]._Type==enSPHNULLPtPair)
				{
					PtiPtr=Region._PtPairList[i]._PtiPtr;
					PtjPtr=Region._PtPairList[i]._PtjPtr;
					KnlPtr=&Region._KnlList[i];

					PtiPtr->_rho+=PtjPtr->_m*KnlPtr->_W;
				}
			}		

			//CSPM-5 修正密度
			for(unsigned int i=0;i<Region._PtList.size();i++)
			{
				if(Region._PtList[i]._Type==enSPHPt&&Region._PtList[i]._Iflag==1&&Region._PtList[i]._CSPMCoef!=0)
				{
					Region._PtList[i]._rho/=Region._PtList[i]._CSPMCoef;
				}
			}
		}
	}


}
