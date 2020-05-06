#include "SPHAVForce.h"

CSPHAVForce::CSPHAVForce()
{
}

CSPHAVForce::~CSPHAVForce()
{
}

void CSPHAVForce::Solve(CRegion &Region)
{
	CSPHPt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;
	double hij;
	double xij2;
	double cij;
	double rhoij;
	double sphi2;
	double bphi;
	double bphi2;
	double uij;
	double vij;
	double wij;
	double uvwxyzij;
	double avforce;
	double avfwx;
	double avfwy;
	double avfwz;
	double avfuvwij;
	double alpha;
	double beta;
	double eta;
	unsigned int IDi,IDj;
  double xij,yij;

	alpha=Region._ControlSPH._AVAlpha;
	beta=Region._ControlSPH._AVBeta;
	eta=Region._ControlSPH._AVEta;

	for (size_t i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair||Region._PtPairList[i]._Type==enSPHNULLPtPair/*||Region._PtPairList[i]._Type==enSPHDumPtPair*/)
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          if(PtiPtr==PtjPtr) continue;

          hij=0.5*(PtiPtr->_h+PtjPtr->_h);
          xij2=Region._PtPairList[i]._driac2;
          //cij=0.5*(PtiPtr->_c+PtjPtr->_c);
          cij=Region._ControlSPH._Cs;
          //cij=10*(0.5*(sqrt(PtiPtr->_u*PtiPtr->_u+PtiPtr->_v*PtiPtr->_v+PtiPtr->_w*PtiPtr->_w)
          //	         +sqrt(PtjPtr->_u*PtjPtr->_u+PtjPtr->_v*PtjPtr->_v+PtjPtr->_w*PtjPtr->_w)));

          rhoij=0.5*(PtiPtr->_rho+PtjPtr->_rho);

          sphi2=eta*eta*hij*hij;
          uij=PtiPtr->_u-PtjPtr->_u;
          vij=PtiPtr->_v-PtjPtr->_v;
          xij=Region._PtPairList[i]._xij;
          yij=Region._PtPairList[i]._yij;
          uvwxyzij=uij*xij+vij*yij;
          bphi=hij*uvwxyzij/(xij2+sphi2);
          bphi2=bphi*bphi;

          if(uvwxyzij<0)
            avforce=(beta*bphi2-alpha*cij*bphi)/rhoij;
          else
            avforce=0;

          avfwx=avforce*KnlPtr->_Wx;
          avfwy=avforce*KnlPtr->_Wy;
          //avfuvwij=avforce*(uij*KnlPtr->_Wx+vij*KnlPtr->_Wy+wij*KnlPtr->_Wz);

          if(Region._ControlSPH._RunMod==1)//显式计算
            {
              PtiPtr->_du-=PtjPtr->_m*avfwx;
              PtiPtr->_dv-=PtjPtr->_m*avfwy;
              //PtiPtr->_de+=0.5*PtjPtr->_m*avfuvwij;

              if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)//PtjPtr->_Type==enSPHPt说明是流体-流体粒子对
                {
                  PtjPtr->_du+=PtiPtr->_m*avfwx;//人工粘度完全对称
                  PtjPtr->_dv+=PtiPtr->_m*avfwy;
                  //PtjPtr->_de+=0.5*PtiPtr->_m*avfuvwij;
                }
            }

          else//隐式计算
            {
              IDi=PtiPtr->_ID2;

              Region._SPHbu[IDi-1]-=PtjPtr->_m*avfwx;
              Region._SPHbv[IDi-1]-=PtjPtr->_m*avfwy;

              if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)//PtjPtr->_Type==enSPHPt说明是流体-流体粒子对
                {
                  IDj=PtjPtr->_ID2;

                  Region._SPHbu[IDj-1]+=PtiPtr->_m*avfwx;
                  Region._SPHbv[IDj-1]+=PtiPtr->_m*avfwy;
                }
            }
        }
    }
}

void CSPHAVForce::Solve2(CRegion &Region)
{
	CSPHPt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;
	double hij;
	double rhoij;
	double rij;
	double xij,yij,zij;
	double uij;
	double vij;
	double wij;
	double Eta;
	unsigned int IDi,IDj;
	double AVkc;//Monaghan1997型人工粘性参数，见SPH-3 P196 式42,k与声速c的乘积
	double PIij;
	double rho,h;
 
	unsigned int LastPtID;

	LastPtID=Region._PtList.size();

	rho=Region._PtList[LastPtID-1]._rho;
	h=Region._PtList[LastPtID-1]._h;

	Eta=Region._PartList[0]._VisK;//暂时用part 0的

	AVkc=(Eta/rho)/(0.134*h);//0.134=15/112
	
	for (size_t i=0;i<Region._PtPairList.size();i++)
	{
		//if(Region._PtPairList[i]._Type==enSPHPtPair||Region._PtPairList[i]._Type==enSPHNULLPtPair||Region._PtPairList[i]._Type==enSPHDumPtPair)
		{
			PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
			PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
			KnlPtr=&Region._KnlList[i];

			if(PtiPtr==PtjPtr) continue;

			hij=0.5*(PtiPtr->_h+PtjPtr->_h);
			rhoij=0.5*(PtiPtr->_rho+PtjPtr->_rho);
			rij=sqrt(Region._PtPairList[i]._driac2);

			xij=Region._PtPairList[i]._xij;
			yij=Region._PtPairList[i]._yij;

			uij=PtiPtr->_u-PtjPtr->_u;
			vij=PtiPtr->_v-PtjPtr->_v;

			PIij=-AVkc*(uij*xij+vij*yij)/(rhoij*rij);

			if(Region._ControlSPH._RunMod==1)//显式计算
			{
				PtiPtr->_du-=PtjPtr->_m*PIij*KnlPtr->_Wx;
				PtiPtr->_dv-=PtjPtr->_m*PIij*KnlPtr->_Wy;

				if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)//PtjPtr->_Type==enSPHPt说明是流体-流体粒子对
				{
					PtjPtr->_du+=PtiPtr->_m*PIij*KnlPtr->_Wx;
					PtjPtr->_dv+=PtiPtr->_m*PIij*KnlPtr->_Wy;
				}
			}

			else//隐式计算
			{
				IDi=PtiPtr->_ID2;

				Region._SPHbu[IDi-1]-=PtjPtr->_m*PIij*KnlPtr->_Wx;
				Region._SPHbv[IDi-1]-=PtjPtr->_m*PIij*KnlPtr->_Wy;

				if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)//PtjPtr->_Type==enSPHPt说明是流体-流体粒子对
				{
					IDj=PtjPtr->_ID2;

					Region._SPHbu[IDj-1]+=PtiPtr->_m*PIij*KnlPtr->_Wx;
					Region._SPHbv[IDj-1]+=PtiPtr->_m*PIij*KnlPtr->_Wy;
				}
			}
		}
	}
}
