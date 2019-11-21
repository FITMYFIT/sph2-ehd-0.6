#include "GetDumProperty.h"

CGetDumProperty::CGetDumProperty()
{

}

CGetDumProperty::~CGetDumProperty()
{

}

void CGetDumProperty::Solve( CRegion &Region )
{
	CBasePt * PtPtr;
	CBasePt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;
	double g;//重力加速度

	double xij,yij;
	double mi,mj;
	double rhoi,rhoj;

	g=Region._ExtForceList[1]->_fy;//有点问题

	//1.将Dummy粒子的压力置0
  // ehd model:set phi=0, 2019.10.15
	for (unsigned int i=0;i!=Region._PtList.size();++i)
	{
		PtPtr=&Region._PtList[i];

		if (PtPtr->_Type==enDumPt||PtPtr->_Type==enEHDBndPt||PtPtr->_Type==enEHDDumPt)
      {
        PtPtr->_p=0.0;

        //PtPtr->_ePhi=0.0;//2019.10.15
      
        //PtPtr->_p=-2*PtPtr->_DistanceBnd*PtPtr->_rho*(g*PtPtr->_Nny);//S Marrone 2011

        PtPtr->_CSPMCoef=0.0;//这里的_CSPMCoef实际上是XYHu　2012 Equ27中的分母（sum(Wwf)）
      }
	}

  //	2.插值得到Dummy粒子的压力值
  //X Y Hu 2012 Equ 27, aw not considered
	for (unsigned int i=0;i!=Region._PtPairList.size();++i)
	{
		if (Region._PtPairList[i]._Type==enSPHDumPtPair
        ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
        ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
      {
        PtiPtr=Region._PtPairList[i]._PtiPtr;
        PtjPtr=Region._PtPairList[i]._PtjPtr;

        KnlPtr=&Region._KnlList[i];

        xij=PtiPtr->_x-PtjPtr->_x;
        yij=PtiPtr->_y-PtjPtr->_y;

        PtjPtr->_CSPMCoef+=KnlPtr->_W;//equ 27 分母

        PtjPtr->_p+=PtiPtr->_p*KnlPtr->_W-g*yij*PtiPtr->_rho*KnlPtr->_W;

        //PtjPtr->_ePhi+=PtiPtr->_ePhi*KnlPtr->_W;
      }
	}

	//double chi=0;//用来求密度的参数，和弱可压缩状态方程的参数一致
	//double gamma=7;
	//double p0=84000;
	//double rho0=1000;
	//for (unsigned int i=0;i!=Region._PtList.size();++i)
	//{
	//	PtPtr=&Region._PtList[i];

	//	if (PtPtr->_Type==enDumPt)
	//	{
	//		if (PtPtr->_CSPMCoef!=0)
	//		{
	//			PtPtr->_p/=PtPtr->_CSPMCoef;

	//			//PtPtr->_rho=rho0*pow((PtPtr->_p-chi)/p0+1,1/gamma);

	//		}
	//	}
	//}




	// //S Marrone 2012Equ 10
  // //calculate the phi in ehd model, for accurate boundary conditions,2019.10.15
	// double cs;//声速
	// cs=Region._ControlSPH._Cs;
	// for (unsigned int i=0;i!=Region._PtPairList.size();++i)
	// {
	// 	if (Region._PtPairList[i]._Type==enSPHDumPtPair)
	// 	{
	// 		PtiPtr=Region._PtPairList[i]._PtiPtr;
	// 		PtjPtr=Region._PtPairList[i]._PtjPtr;

	// 		KnlPtr=&Region._KnlList[i];

	// 		yij=PtiPtr->_y-PtjPtr->_y;

	// 		PtjPtr->_CSPMCoef+=KnlPtr->_W;

	// 		PtjPtr->_p+=(PtiPtr->_p-PtiPtr->_rho*cs*(PtiPtr->_u*PtjPtr->_Nnx+PtiPtr->_v*PtjPtr->_Nny)-PtiPtr->_rho*g*yij)*KnlPtr->_W;
	// 	}
	// }

	for (unsigned int i=0;i!=Region._PtList.size();++i)
	{
		PtPtr=&Region._PtList[i];

		if (PtPtr->_Type==enDumPt&&PtPtr->_CSPMCoef!=0)
		{
			PtPtr->_p/=PtPtr->_CSPMCoef;

      PtPtr->_ePhi/=PtPtr->_CSPMCoef;
		}
	}



	//非滑移边界条件，插值得到dummy粒子的速度
	//将dummy粒子的速度置0
	for (unsigned int i=0;i!=Region._PtList.size();++i)
	{
		PtPtr=&Region._PtList[i];

		if (PtPtr->_Type==enDumPt/*&&PtPtr->_PID==1*/)
		{
			PtPtr->_u=0.0;
			PtPtr->_v=0.0;

			PtPtr->_CSPMCoef=0.0;
		}
	}

	//插值得到dummy粒子的速度，XY HU2012 Equ22
	for (unsigned int i=0;i!=Region._PtPairList.size();++i)
	{
		if (Region._PtPairList[i]._Type==enSPHDumPtPair)
		{
			PtiPtr=Region._PtPairList[i]._PtiPtr;
			PtjPtr=Region._PtPairList[i]._PtjPtr;

			KnlPtr=&Region._KnlList[i];

			//if(PtjPtr->_PID==1)
			{
				PtjPtr->_CSPMCoef+=KnlPtr->_W;//equ 22 分母

				PtjPtr->_u+=PtiPtr->_u*KnlPtr->_W;
				PtjPtr->_v+=PtiPtr->_v*KnlPtr->_W;
			}
		}
	}

	for (unsigned int i=0;i!=Region._PtList.size();++i)
	{
		PtPtr=&Region._PtList[i];

		if (PtPtr->_Type==enDumPt||PtPtr->_Type==enEHDDumPt||PtPtr->_Type==enEHDBndPt)
      {
        if (PtPtr->_CSPMCoef!=0)
          {
            PtPtr->_u/=PtPtr->_CSPMCoef;
            PtPtr->_v/=PtPtr->_CSPMCoef;

            PtPtr->_u=2*PtPtr->_uwall-PtPtr->_u;
            PtPtr->_v=2*PtPtr->_vwall-PtPtr->_v;

            PtPtr->_CSPMCoef=0.0;
          }
      }
	}

}
