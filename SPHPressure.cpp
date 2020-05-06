#include "SPHPressure.h"

CSPHPressure::CSPHPressure()
{
}

CSPHPressure::~CSPHPressure()
{
}

void CSPHPressure::Solve(CRegion &Region)
{
	CSPHPt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;

	double norm;
	double normx,normy;
	double uij,vij,wij;
	unsigned int IDi,IDj;
	double wxc,wyc,wzc;//修正核函数梯度
	double Vi,Vj;//i j 体积
	double paverg;

	for (size_t i=0;i<Region._PtPairList.size();i++)
	{
		if(Region._PtPairList[i]._Type==enSPHPtPair)
      {
        PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
        PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
        KnlPtr=&Region._KnlList[i];

        //norm=PtiPtr->_p/(PtiPtr->_rho*PtiPtr->_rho)+PtjPtr->_p/(PtjPtr->_rho*PtjPtr->_rho);
        norm=(PtiPtr->_p+PtjPtr->_p)/(PtiPtr->_rho*PtjPtr->_rho);

        ////X Y Hu 2012 Equ(7)
        //Vi=PtiPtr->_m/PtiPtr->_rho;
        //Vj=PtjPtr->_m/PtjPtr->_rho;
        //paverg=(PtjPtr->_rho*PtiPtr->_p+PtiPtr->_rho*PtjPtr->_p)/(PtiPtr->_rho+PtjPtr->_rho);//X Y Hu 2012 Equ(8)

        //normx=(Vi*Vi+Vj*Vj)*paverg*KnlPtr->_Wx;
        //normy=(Vi*Vi+Vj*Vj)*paverg*KnlPtr->_Wy;

        uij=PtiPtr->_u-PtjPtr->_u;
        vij=PtiPtr->_v-PtjPtr->_v;

        if(Region._ControlSPH._RunMod==1)//显式计算
          {
            //PtiPtr->_du-=1/PtiPtr->_m*normx;
            //PtiPtr->_dv-=1/PtiPtr->_m*normy;

            PtiPtr->_du-=PtjPtr->_m*norm*KnlPtr->_Wx;
            PtiPtr->_dv-=PtjPtr->_m*norm*KnlPtr->_Wy;

            if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)
              {
                //PtjPtr->_du+=1/PtjPtr->_m*normx;
                //PtjPtr->_dv+=1/PtjPtr->_m*normy;

                PtjPtr->_du+=PtiPtr->_m*norm*KnlPtr->_Wx;
                PtjPtr->_dv+=PtiPtr->_m*norm*KnlPtr->_Wy;
              }
          }

        else//隐式计算
          {
            IDi=PtiPtr->_ID2;
					
            //_GetKnlListP.GetCrctKnlGrad(Region,PtiPtr,KnlPtr,wxc,wyc,wzc);//求出修正核函数

            //Region._SPHbu[IDi-1]-=PtjPtr->_m*norm*wxc;
            //Region._SPHbv[IDi-1]-=PtjPtr->_m*norm*wyc;

            Region._SPHbu[IDi-1]-=PtjPtr->_m*norm*KnlPtr->_Wx;
            Region._SPHbv[IDi-1]-=PtjPtr->_m*norm*KnlPtr->_Wy;

            ////粒子光滑长度不一致时的求解公式
            //normx=PtjPtr->_m/(PtiPtr->_rho*PtjPtr->_rho)*(PtiPtr->_p*KnlPtr->_Wxj-PtjPtr->_p*KnlPtr->_Wx);//R Vacondio 2013 Equ.6
            //normy=PtjPtr->_m/(PtiPtr->_rho*PtjPtr->_rho)*(PtiPtr->_p*KnlPtr->_Wyj-PtjPtr->_p*KnlPtr->_Wy);

            //Region._SPHbu[IDi-1]-=normx;
            //Region._SPHbv[IDi-1]-=normy;

            if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)
              {
                IDj=PtjPtr->_ID2;
				
                //_GetKnlListP.GetCrctKnlGrad(Region,PtjPtr,KnlPtr,wxc,wyc,wzc);//求出修正核函数

                //Region._SPHbu[IDj-1]+=PtiPtr->_m*norm*wxc;
                //Region._SPHbv[IDj-1]+=PtiPtr->_m*norm*wyc;

                Region._SPHbu[IDj-1]+=PtiPtr->_m*norm*KnlPtr->_Wx;
                Region._SPHbv[IDj-1]+=PtiPtr->_m*norm*KnlPtr->_Wy;

                //Region._SPHbu[IDj-1]+=normx;
                //Region._SPHbv[IDj-1]+=normy;
              }
          }
		}
	}
}
