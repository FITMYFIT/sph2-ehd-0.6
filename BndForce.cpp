#include "BndForce.h"

CBndForce::CBndForce()
{
}

CBndForce::~CBndForce()
{
}

void CBndForce::Solve(CRegion &Region)
{
	CBasePt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;
	double uij,vij;
	double njx,njy;
	double xij,yij;
	double bndnorm;
	double cs;
	double hij;
	unsigned int IDi,IDj;

	double pi,pj;
	double mi,mj;
	double rhoi,rhoj;
	double modulev;

	//for (size_t i=0;i<Region._PtPairList.size();i++)
	//{
	//	if (Region._PtPairList[i]._Type==enSPHBndPtPair/*||Region._PtPairList[i]._Type==enSPHNULLPtPair*/)
	//	{
	//		PtiPtr=Region._PtPairList[i]._PtiPtr;
	//		PtjPtr=Region._PtPairList[i]._PtjPtr;

	//		KnlPtr=&Region._KnlList[i];
	//		
	//		IDi=PtiPtr->_ID2;
	//		IDj=PtjPtr->_ID2;

	//		xij=PtiPtr->_x-PtjPtr->_x;
	//		yij=PtiPtr->_y-PtjPtr->_y;

	//		//2.边界力项
	//		uij=PtiPtr->_u-PtjPtr->_u;
	//		vij=PtiPtr->_v-PtjPtr->_v;

	//		njx=PtjPtr->_Nnx;
	//		njy=PtjPtr->_Nny;

	//		if((uij*njx+vij*njy)<0)//韩边界施加条件
	//		{
	//			cs=Region._ControlSPH._Cs;
	//			hij=(PtiPtr->_h+PtjPtr->_h)/2;
	//			//hij=PtiPtr->_h;

	//			bndnorm=-0.01*cs*cs*_OperatorB.LMin(uij*njx+vij*njy,-1)*hij*hij*KnlPtr->_W/fabs(xij*njx+yij*njy);


	//			if (Region._ControlSPH._RunMod==1)//显式算法
	//			{
	//				PtiPtr->_du+=bndnorm*njx;
	//				PtiPtr->_dv+=bndnorm*njy;
	//			}

	//			else
	//			{
	//				Region._SPHbu[IDi-1]+=bndnorm*njx;
	//				Region._SPHbv[IDi-1]+=bndnorm*njy;
	//			}
	//		}

	//	}
	//}


  //	X Y Hu压力边界施加方法
	for (unsigned int i=0;i!=Region._PtPairList.size();++i)
	{
		if (Region._PtPairList[i]._Type==enSPHDumPtPair
        ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
        ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
      {		
        PtiPtr=Region._PtPairList[i]._PtiPtr;
        PtjPtr=Region._PtPairList[i]._PtjPtr;

        KnlPtr=&Region._KnlList[i];

        njx=PtjPtr->_Nnx;
        njy=PtjPtr->_Nny;

        pi=PtiPtr->_p;
        pj=PtjPtr->_p;

        mi=PtiPtr->_m;
        mj=PtjPtr->_m;

        rhoi=PtiPtr->_rho;
        rhoj=PtjPtr->_rho;

        bndnorm=mj*(pi+pj)/rhoi/rhoj;

        PtiPtr->_du-=bndnorm*KnlPtr->_Wx;
        PtiPtr->_dv-=bndnorm*KnlPtr->_Wy;

        PtiPtr->_AccBndx-=bndnorm*KnlPtr->_Wx;
        PtiPtr->_AccBndy-=bndnorm*KnlPtr->_Wy;
      }		
	}


	// //新压力边界施加方法
	// //fi=sum(mj/rhoj*(Pi+Pj)/rhoi/rhoj*(Wx*Nnx+Wy*Nny))
	// for (unsigned int i=0;i!=Region._PtPairList.size();++i)
	// {
	// 	if (Region._PtPairList[i]._Type==enSPHDumPtPair
  //       ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
  //       ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
  //     {		
  //       PtiPtr=Region._PtPairList[i]._PtiPtr;
  //       PtjPtr=Region._PtPairList[i]._PtjPtr;

  //       KnlPtr=&Region._KnlList[i];

  //       njx=PtjPtr->_Nnx;
  //       njy=PtjPtr->_Nny;

  //       pi=PtiPtr->_p;
  //       pj=PtjPtr->_p;

  //       mi=PtiPtr->_m;
  //       mj=PtjPtr->_m;

  //       rhoi=PtiPtr->_rho;
  //       rhoj=PtjPtr->_rho;

  //       uij=PtiPtr->_u-PtjPtr->_u;
  //       vij=PtiPtr->_v-PtjPtr->_v;

  //       modulev=sqrt(PtiPtr->_u*PtiPtr->_u+PtiPtr->_v*PtiPtr->_v);

  //       //if (fabs(uij*PtjPtr->_Nnx+vij*PtjPtr->_Nny)>=0.01*modulev)
  //       {
  //         if (pi+pj>0)
  //           {
  //             bndnorm=mj*(pi+pj)/rhoi/rhoj*(KnlPtr->_Wx*njx+KnlPtr->_Wy*njy);
  //             //bndnorm=mj*(pi/(rhoi*rhoi)+pj/(rhoj*rhoj))*(KnlPtr->_Wx*njx+KnlPtr->_Wy*njy);

  //             PtiPtr->_du-=bndnorm*njx;
  //             PtiPtr->_dv-=bndnorm*njy;

  //             PtiPtr->_AccBndx-=bndnorm*njx;
  //             PtiPtr->_AccBndy-=bndnorm*njy;
  //           }
  //       }			
  //     }		
  //}

}
