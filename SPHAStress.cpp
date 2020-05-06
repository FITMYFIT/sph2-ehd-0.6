#include "SPHAStress.h"
#include <cmath>

CSPHAStress::CSPHAStress()
{
}

CSPHAStress::~CSPHAStress()
{
}

void CSPHAStress::Solve(CRegion &Region)
{
	CSPHPt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;
	double fij;
	double n;
	double rixx;
	double riyy;
	double rizz;
	double rjxx;
	double rjyy;
	double rjzz;
	double normx;
	double normy;
	double normz;
	double normwx;
	double normwy;
	double normwz;
	double epsilon1;
	double epsilon2;
	double rhoi2;
	double rhoj2;
	double w0;

	unsigned int IDi,IDj;

	epsilon1=Region._ControlSPH._ASEpsilon1;
	epsilon2=Region._ControlSPH._ASEpsilon2;

	for (size_t i=0;i<Region._PtPairList.size();i++)
	{
		if(Region._PtPairList[i]._Type!=enSPHBndPtPair/*&&Region._PtPairList[i]._Type!=enSPHDumPtPair*/)//除了流体-边界粒子对外，均计算人工应力
		{
			PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
			PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;

			if(PtiPtr==PtjPtr) continue;

			KnlPtr=&Region._KnlList[i];
			_GetKnlList.GetW0(Region._ControlSPH._ASDeltaD,0.5*(PtiPtr->_h+PtjPtr->_h),w0);
			fij=KnlPtr->_W/w0;
			n=4;

			rhoi2=PtiPtr->_rho*PtiPtr->_rho;
			rhoj2=PtjPtr->_rho*PtjPtr->_rho;

			if(PtiPtr->_p<0)
			{
				rixx=-epsilon1*(PtiPtr->_p-PtiPtr->_sxx)/rhoi2;
				riyy=-epsilon1*(PtiPtr->_p-PtiPtr->_syy)/rhoi2;
			}
			else
			{
				rixx=epsilon2*(PtiPtr->_p-PtiPtr->_sxx)/rhoi2;
				riyy=epsilon2*(PtiPtr->_p-PtiPtr->_syy)/rhoi2;
			}

			if(PtjPtr->_p<0)
			{
				rjxx=-epsilon1*(PtjPtr->_p-PtjPtr->_sxx)/rhoj2;
				rjyy=-epsilon1*(PtjPtr->_p-PtjPtr->_syy)/rhoj2;
			}
			else
			{
				rjxx=epsilon2*(PtjPtr->_p-PtjPtr->_sxx)/rhoj2;
				rjyy=epsilon2*(PtjPtr->_p-PtjPtr->_syy)/rhoj2;
			}

			normx=pow(fij,n)*(rixx+rjxx);
			normy=pow(fij,n)*(riyy+rjyy);

			normwx=normx*KnlPtr->_Wx;
			normwy=normy*KnlPtr->_Wy;

			if(Region._ControlSPH._RunMod==1)//显式计算
			{
				PtiPtr->_du-=PtjPtr->_m*normwx;
				PtiPtr->_dv-=PtjPtr->_m*normwy;

				if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)//PtjPtr->_Type==enSPHPt说明是流体-流体粒子对
				{
					PtjPtr->_du+=PtiPtr->_m*normwx;
					PtjPtr->_dv+=PtiPtr->_m*normwy;
				}
			}

			else//隐式计算
			{
				IDi=PtiPtr->_ID2;

				Region._SPHbu[IDi-1]-=PtjPtr->_m*normwx;
				Region._SPHbv[IDi-1]-=PtjPtr->_m*normwy;

				if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)//PtjPtr->_Type==enSPHPt说明是流体-流体粒子对
				{
					IDj=PtjPtr->_ID2;

					Region._SPHbu[IDj-1]+=PtiPtr->_m*normwx;
					Region._SPHbv[IDj-1]+=PtiPtr->_m*normwy;
				}
			}
		}
	}
}
