#include "SPHShearingForce.h"

CSPHShearingForce::CSPHShearingForce()
{
}

CSPHShearingForce::~CSPHShearingForce()
{
}

void CSPHShearingForce::UpdateDs(CRegion &Region)
{
//	CSPHPt * PtiPtr,* PtjPtr;
//	CKnl * KnlPtr;
//
//	CPart * PartPtr;
//	CSPHPt * SPHPtPtr;
//
//	double G;
//	double norm;
//	double uij,vij,wij;
//
//	G=Region._ControlSPH._G;
//
//	for (size_t i=0;i<Region._PtPairList.size();i++)
//	{
//		if(Region._PtPairList[i]._Type==enSPHPtPair)
//		{
//			PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
//			PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
//			KnlPtr=&Region._KnlList[i];
//
//			uij=PtiPtr->_u-PtjPtr->_u;
//			vij=PtiPtr->_v-PtjPtr->_v;
//			wij=PtiPtr->_w-PtjPtr->_w;
//
//			PtiPtr->_depsilonxx-=(PtjPtr->_m/PtjPtr->_rho)*(uij*KnlPtr->_Wx);
//			PtiPtr->_depsilonyy-=(PtjPtr->_m/PtjPtr->_rho)*(vij*KnlPtr->_Wy);
//			PtiPtr->_depsilonzz-=(PtjPtr->_m/PtjPtr->_rho)*(wij*KnlPtr->_Wz);
//			PtiPtr->_depsilonxy-=0.5*(PtjPtr->_m/PtjPtr->_rho)*(uij*KnlPtr->_Wy+vij*KnlPtr->_Wx);
//			PtiPtr->_depsilonxz-=0.5*(PtjPtr->_m/PtjPtr->_rho)*(uij*KnlPtr->_Wz+wij*KnlPtr->_Wx);
//			PtiPtr->_depsilonyz-=0.5*(PtjPtr->_m/PtjPtr->_rho)*(vij*KnlPtr->_Wz+wij*KnlPtr->_Wy);
//
//			PtiPtr->_Rxy-=0.5*(PtjPtr->_m/PtjPtr->_rho)*(uij*KnlPtr->_Wy-vij*KnlPtr->_Wx);
//			PtiPtr->_Rxz-=0.5*(PtjPtr->_m/PtjPtr->_rho)*(uij*KnlPtr->_Wz-wij*KnlPtr->_Wx);
//			PtiPtr->_Ryz-=0.5*(PtjPtr->_m/PtjPtr->_rho)*(vij*KnlPtr->_Wz-wij*KnlPtr->_Wy);
//			if(PtiPtr!=PtjPtr)
//			{
//				PtjPtr->_depsilonxx-=(PtiPtr->_m/PtiPtr->_rho)*(uij*KnlPtr->_Wx);
//				PtjPtr->_depsilonyy-=(PtiPtr->_m/PtiPtr->_rho)*(vij*KnlPtr->_Wy);
//				PtjPtr->_depsilonzz-=(PtiPtr->_m/PtiPtr->_rho)*(wij*KnlPtr->_Wz);
//				PtjPtr->_depsilonxy-=0.5*(PtiPtr->_m/PtiPtr->_rho)*(uij*KnlPtr->_Wy+vij*KnlPtr->_Wx);
//				PtjPtr->_depsilonxz-=0.5*(PtiPtr->_m/PtiPtr->_rho)*(uij*KnlPtr->_Wz+wij*KnlPtr->_Wx);
//				PtjPtr->_depsilonyz-=0.5*(PtiPtr->_m/PtiPtr->_rho)*(vij*KnlPtr->_Wz+wij*KnlPtr->_Wy);
//
//				PtjPtr->_Rxy-=0.5*(PtiPtr->_m/PtiPtr->_rho)*(uij*KnlPtr->_Wy-vij*KnlPtr->_Wx);
//				PtjPtr->_Rxz-=0.5*(PtiPtr->_m/PtiPtr->_rho)*(uij*KnlPtr->_Wz-wij*KnlPtr->_Wx);
//				PtjPtr->_Ryz-=0.5*(PtiPtr->_m/PtiPtr->_rho)*(vij*KnlPtr->_Wz-wij*KnlPtr->_Wy);
//			}
//		}
//	}
//
//	for (size_t i=0;i<Region._PartList.size();i++)
//	{
//		PartPtr=&Region._PartList[i];
//		if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
//			||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
//		{
//			for(size_t j=0;j<PartPtr->_PartPtList.size();j++)
//			{
//				SPHPtPtr=PartPtr->_PartPtList[j];
//				norm=0.666667*(SPHPtPtr->_depsilonxx+SPHPtPtr->_depsilonyy+SPHPtPtr->_depsilonzz);
//				SPHPtPtr->_dsxx=G*(2*SPHPtPtr->_depsilonxx-norm)+2*SPHPtPtr->_sxy*SPHPtPtr->_Rxy+2*SPHPtPtr->_sxz*SPHPtPtr->_Rxz;
//				SPHPtPtr->_dsyy=G*(2*SPHPtPtr->_depsilonyy-norm)-2*SPHPtPtr->_sxy*SPHPtPtr->_Rxy+2*SPHPtPtr->_syz*SPHPtPtr->_Ryz;
//				SPHPtPtr->_dszz=G*(2*SPHPtPtr->_depsilonzz-norm)-2*SPHPtPtr->_sxz*SPHPtPtr->_Rxz-2*SPHPtPtr->_syz*SPHPtPtr->_Ryz;
//				SPHPtPtr->_dsxy=2*G*SPHPtPtr->_depsilonxy+SPHPtPtr->_Rxy*(SPHPtPtr->_syy-SPHPtPtr->_sxx)+SPHPtPtr->_Rxz*SPHPtPtr->_syz+SPHPtPtr->_Ryz*SPHPtPtr->_sxz;
//				SPHPtPtr->_dsxz=2*G*SPHPtPtr->_depsilonxz+SPHPtPtr->_Rxz*(SPHPtPtr->_szz-SPHPtPtr->_sxx)+SPHPtPtr->_Rxy*SPHPtPtr->_syz-SPHPtPtr->_Ryz*SPHPtPtr->_sxy;
//				SPHPtPtr->_dsyz=2*G*SPHPtPtr->_depsilonyz+SPHPtPtr->_Ryz*(SPHPtPtr->_szz-SPHPtPtr->_syy)-SPHPtPtr->_Rxz*SPHPtPtr->_sxy-SPHPtPtr->_Rxy*SPHPtPtr->_sxz;
//			}
//		}
//	}
//}
//
//void CSPHShearingForce::UpdateS(CRegion &Region,double DeltaT)
//{
//	CPart * PartPtr;
//	CSPHPt * SPHPtPtr;
//
//	double YieldStress;
//
//	YieldStress=Region._ControlSPH._YieldStress;
//
//	for (size_t i=0;i<Region._PartList.size();i++)
//	{
//		PartPtr=&Region._PartList[i];
//		if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
//			||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
//		{
//			for(size_t j=0;j<PartPtr->_PartPtList.size();j++)
//			{
//				SPHPtPtr=PartPtr->_PartPtList[j];
//				SPHPtPtr->_sxx+=SPHPtPtr->_dsxx*DeltaT;
//				SPHPtPtr->_sxy+=SPHPtPtr->_dsxy*DeltaT;
//				SPHPtPtr->_sxz+=SPHPtPtr->_dsxz*DeltaT;
//				SPHPtPtr->_syy+=SPHPtPtr->_dsyy*DeltaT;
//				SPHPtPtr->_syz+=SPHPtPtr->_dsyz*DeltaT;
//				SPHPtPtr->_szz+=SPHPtPtr->_dszz*DeltaT;
//
//				SPHPtPtr->_J2=sqrt(1.5*(SPHPtPtr->_sxx*SPHPtPtr->_sxx+SPHPtPtr->_syy*SPHPtPtr->_syy+SPHPtPtr->_szz*SPHPtPtr->_szz
//					+2*SPHPtPtr->_sxy*SPHPtPtr->_sxy+2*SPHPtPtr->_sxz*SPHPtPtr->_sxz+2*SPHPtPtr->_syz*SPHPtPtr->_syz));
//				if(SPHPtPtr->_J2>YieldStress)
//				{
//					double norm;
//					norm=YieldStress/SPHPtPtr->_J2;
//					SPHPtPtr->_sxx*=norm;
//					SPHPtPtr->_sxy*=norm;
//					SPHPtPtr->_sxz*=norm;
//					SPHPtPtr->_syy*=norm;
//					SPHPtPtr->_syz*=norm;
//					SPHPtPtr->_szz*=norm;
//				}
//			}
//		}
//	}
//}
//
//void CSPHShearingForce::Solve(CRegion &Region,double DeltaT)
//{
//	UpdateDs(Region);
//	UpdateS(Region,DeltaT);
//
//	CSPHPt * PtiPtr,* PtjPtr;
//	CKnl * KnlPtr;
//
//	double normx;
//	double normy;
//	double normz;
//
//
//	for (size_t i=0;i<Region._PtPairList.size();i++)
//	{
//		if(Region._PtPairList[i]._Type==enSPHPtPair)
//		{
//			PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
//			PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
//			KnlPtr=&Region._KnlList[i];
//
//			normx=((PtiPtr->_sxx+PtjPtr->_sxx)*KnlPtr->_Wx+(PtiPtr->_sxy+PtjPtr->_sxy)*KnlPtr->_Wy+(PtiPtr->_sxz+PtjPtr->_sxz)*KnlPtr->_Wz)/(PtiPtr->_rho*PtjPtr->_rho);
//			normy=((PtiPtr->_sxy+PtjPtr->_sxy)*KnlPtr->_Wx+(PtiPtr->_syy+PtjPtr->_syy)*KnlPtr->_Wy+(PtiPtr->_syz+PtjPtr->_syz)*KnlPtr->_Wz)/(PtiPtr->_rho*PtjPtr->_rho);
//			normz=((PtiPtr->_sxz+PtjPtr->_sxz)*KnlPtr->_Wx+(PtiPtr->_syz+PtjPtr->_syz)*KnlPtr->_Wy+(PtiPtr->_szz+PtjPtr->_szz)*KnlPtr->_Wz)/(PtiPtr->_rho*PtjPtr->_rho);
//
//			PtiPtr->_du+=PtjPtr->_m*normx;
//			PtiPtr->_dv+=PtjPtr->_m*normy;
//			PtiPtr->_dw+=PtjPtr->_m*normz;
//			//PtiPtr->_de+=???
//			if(PtiPtr!=PtjPtr)
//			{
//				PtjPtr->_du-=PtiPtr->_m*normx;
//				PtjPtr->_dv-=PtiPtr->_m*normy;
//				PtjPtr->_dw-=PtiPtr->_m*normz;
//				//PtiPtr->_de+=???
//			}
//		}
//	}
}
