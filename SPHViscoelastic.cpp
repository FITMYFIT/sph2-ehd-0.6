#include "SPHViscoelastic.h"

CSPHViscoelastic::CSPHViscoelastic()
{
}

CSPHViscoelastic::~CSPHViscoelastic()
{
}

void CSPHViscoelastic::UpdateDs(CRegion &Region)
{
	//CSPHPt * PtiPtr,* PtjPtr;
	//CKnl * KnlPtr;

	//CPart * PartPtr;
	//CSPHPt * SPHPtPtr;

	//double visetap;
	//double vislambda;
	//double visg;
	//double uij,vij,wij;

	//visetap=Region._ControlSPH._VisEtap;
	//visg=Region._ControlSPH._Visg;
	//vislambda=Region._ControlSPH._VisLambda;

	//for (size_t i=0;i<Region._PtPairList.size();i++)
	//{
	//	if(Region._PtPairList[i]._Type==enSPHPtPair)
	//	{
	//		PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
	//		PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
	//		KnlPtr=&Region._KnlList[i];

	//		uij=PtiPtr->_u-PtjPtr->_u;
	//		vij=PtiPtr->_v-PtjPtr->_v;
	//		wij=PtiPtr->_w-PtjPtr->_w;

	//		PtiPtr->_kxx-=(PtjPtr->_m/PtjPtr->_rho)*uij*KnlPtr->_Wx;
	//		PtiPtr->_kxy-=(PtjPtr->_m/PtjPtr->_rho)*uij*KnlPtr->_Wy;
	//		PtiPtr->_kxz-=(PtjPtr->_m/PtjPtr->_rho)*uij*KnlPtr->_Wz;
	//		PtiPtr->_kyx-=(PtjPtr->_m/PtjPtr->_rho)*vij*KnlPtr->_Wx;
	//		PtiPtr->_kyy-=(PtjPtr->_m/PtjPtr->_rho)*vij*KnlPtr->_Wy;
	//		PtiPtr->_kyz-=(PtjPtr->_m/PtjPtr->_rho)*vij*KnlPtr->_Wz;
	//		PtiPtr->_kzx-=(PtjPtr->_m/PtjPtr->_rho)*wij*KnlPtr->_Wx;
	//		PtiPtr->_kzy-=(PtjPtr->_m/PtjPtr->_rho)*wij*KnlPtr->_Wy;
	//		PtiPtr->_kzz-=(PtjPtr->_m/PtjPtr->_rho)*wij*KnlPtr->_Wz;

	//		if(PtiPtr!=PtjPtr)
	//		{
	//			PtjPtr->_kxx-=(PtiPtr->_m/PtiPtr->_rho)*uij*KnlPtr->_Wx;
	//			PtjPtr->_kxy-=(PtiPtr->_m/PtiPtr->_rho)*uij*KnlPtr->_Wy;
	//			PtjPtr->_kxz-=(PtiPtr->_m/PtiPtr->_rho)*uij*KnlPtr->_Wz;
	//			PtjPtr->_kyx-=(PtiPtr->_m/PtiPtr->_rho)*vij*KnlPtr->_Wx;
	//			PtjPtr->_kyy-=(PtiPtr->_m/PtiPtr->_rho)*vij*KnlPtr->_Wy;
	//			PtjPtr->_kyz-=(PtiPtr->_m/PtiPtr->_rho)*vij*KnlPtr->_Wz;
	//			PtjPtr->_kzx-=(PtiPtr->_m/PtiPtr->_rho)*wij*KnlPtr->_Wx;
	//			PtjPtr->_kzy-=(PtiPtr->_m/PtiPtr->_rho)*wij*KnlPtr->_Wy;
	//			PtjPtr->_kzz-=(PtiPtr->_m/PtiPtr->_rho)*wij*KnlPtr->_Wz;
	//		}
	//	}
	//}

	//for (size_t i=0;i<Region._PartList.size();i++)
	//{
	//	PartPtr=&Region._PartList[i];
	//	if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
	//		||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
	//	{
	//		for(size_t j=0;j<PartPtr->_PartPtList.size();j++)
	//		{
	//			SPHPtPtr=PartPtr->_PartPtList[j];

	//			SPHPtPtr->_depsilonxx=SPHPtPtr->_kxx;
	//			SPHPtPtr->_depsilonyy=SPHPtPtr->_kyy;
	//			SPHPtPtr->_depsilonzz=SPHPtPtr->_kzz;
	//			SPHPtPtr->_depsilonxy=0.5*(SPHPtPtr->_kxy+SPHPtPtr->_kyx);
	//			SPHPtPtr->_depsilonxz=0.5*(SPHPtPtr->_kxz+SPHPtPtr->_kzx);
	//			SPHPtPtr->_depsilonyz=0.5*(SPHPtPtr->_kyz+SPHPtPtr->_kzy);

	//			if(vislambda==0.0)
	//			{
	//				SPHPtPtr->_dsxx=0.0;
	//				SPHPtPtr->_dsyy=0.0;
	//				SPHPtPtr->_dszz=0.0;
	//				SPHPtPtr->_dsxy=0.0;
	//				SPHPtPtr->_dsxz=0.0;
	//				SPHPtPtr->_dsyz=0.0;

	//				SPHPtPtr->_sxx=2*visetap*SPHPtPtr->_depsilonxx;
	//				SPHPtPtr->_sxy=2*visetap*SPHPtPtr->_depsilonxy;
	//				SPHPtPtr->_sxz=2*visetap*SPHPtPtr->_depsilonxz;
	//				SPHPtPtr->_syy=2*visetap*SPHPtPtr->_depsilonyy;
	//				SPHPtPtr->_syz=2*visetap*SPHPtPtr->_depsilonyz;
	//				SPHPtPtr->_szz=2*visetap*SPHPtPtr->_depsilonzz;
	//			}
	//			else
	//			{
	//				SPHPtPtr->_dsxx=2*SPHPtPtr->_kxx*SPHPtPtr->_sxx+2*SPHPtPtr->_kxy*SPHPtPtr->_sxy+2*SPHPtPtr->_kxz*SPHPtPtr->_sxz
	//					-(visg/vislambda)*SPHPtPtr->_sxx+2*(visetap/vislambda)*SPHPtPtr->_kxx;

	//				SPHPtPtr->_dsyy=2*SPHPtPtr->_kyx*SPHPtPtr->_sxy+2*SPHPtPtr->_kyy*SPHPtPtr->_syy+2*SPHPtPtr->_kyz*SPHPtPtr->_syz
	//					-(visg/vislambda)*SPHPtPtr->_syy+2*(visetap/vislambda)*SPHPtPtr->_kyy;

	//				SPHPtPtr->_dszz=2*SPHPtPtr->_kzx*SPHPtPtr->_sxz+2*SPHPtPtr->_kzy*SPHPtPtr->_syz+2*SPHPtPtr->_kzz*SPHPtPtr->_szz
	//					-(visg/vislambda)*SPHPtPtr->_szz+2*(visetap/vislambda)*SPHPtPtr->_kzz;

	//				SPHPtPtr->_dsxy=SPHPtPtr->_kxx*SPHPtPtr->_sxy+SPHPtPtr->_kyx*SPHPtPtr->_sxx+SPHPtPtr->_kxy*SPHPtPtr->_syy
	//					+SPHPtPtr->_kyy*SPHPtPtr->_sxy+SPHPtPtr->_kxz*SPHPtPtr->_syz+SPHPtPtr->_kyz*SPHPtPtr->_sxz
	//					-(visg/vislambda)*SPHPtPtr->_sxy+2*(visetap/vislambda)*SPHPtPtr->_depsilonxy;

	//				SPHPtPtr->_dsxz=SPHPtPtr->_kxx*SPHPtPtr->_sxz+SPHPtPtr->_kzx*SPHPtPtr->_sxx+SPHPtPtr->_kxy*SPHPtPtr->_syz
	//					+SPHPtPtr->_kzy*SPHPtPtr->_sxy+SPHPtPtr->_kxz*SPHPtPtr->_szz+SPHPtPtr->_kzz*SPHPtPtr->_sxz
	//					-(visg/vislambda)*SPHPtPtr->_sxz+2*(visetap/vislambda)*SPHPtPtr->_depsilonxz;

	//				SPHPtPtr->_dsyz=SPHPtPtr->_kyx*SPHPtPtr->_sxz+SPHPtPtr->_kzx*SPHPtPtr->_sxy+SPHPtPtr->_kyy*SPHPtPtr->_syz
	//					+SPHPtPtr->_kzy*SPHPtPtr->_syy+SPHPtPtr->_kyz*SPHPtPtr->_szz+SPHPtPtr->_kzz*SPHPtPtr->_syz
	//					-(visg/vislambda)*SPHPtPtr->_syz+2*(visetap/vislambda)*SPHPtPtr->_depsilonyz;
	//			}
	//		}
	//	}
	//}
}

void CSPHViscoelastic::UpdateS(CRegion &Region,double DeltaT)
{
	//CPart * PartPtr;
	//CSPHPt * SPHPtPtr;

	//double visetas;
	//visetas=Region._ControlSPH._VisEtas;

	//for (size_t i=0;i<Region._PartList.size();i++)
	//{
	//	PartPtr=&Region._PartList[i];
	//	if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
	//		||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
	//	{
	//		for(size_t j=0;j<PartPtr->_PartPtList.size();j++)
	//		{
	//			SPHPtPtr=PartPtr->_PartPtList[j];
	//			SPHPtPtr->_sxx+=SPHPtPtr->_dsxx*DeltaT;
	//			SPHPtPtr->_sxy+=SPHPtPtr->_dsxy*DeltaT;
	//			SPHPtPtr->_sxz+=SPHPtPtr->_dsxz*DeltaT;
	//			SPHPtPtr->_syy+=SPHPtPtr->_dsyy*DeltaT;
	//			SPHPtPtr->_syz+=SPHPtPtr->_dsyz*DeltaT;
	//			SPHPtPtr->_szz+=SPHPtPtr->_dszz*DeltaT;

	//			SPHPtPtr->_sigmaxx=2*visetas*SPHPtPtr->_depsilonxx+SPHPtPtr->_sxx;
	//			SPHPtPtr->_sigmaxy=2*visetas*SPHPtPtr->_depsilonxy+SPHPtPtr->_sxy;
	//			SPHPtPtr->_sigmaxz=2*visetas*SPHPtPtr->_depsilonxz+SPHPtPtr->_sxz;
	//			SPHPtPtr->_sigmayy=2*visetas*SPHPtPtr->_depsilonyy+SPHPtPtr->_syy;
	//			SPHPtPtr->_sigmayz=2*visetas*SPHPtPtr->_depsilonyz+SPHPtPtr->_syz;
	//			SPHPtPtr->_sigmazz=2*visetas*SPHPtPtr->_depsilonzz+SPHPtPtr->_szz;
	//		}
	//	}
	//}
}

void CSPHViscoelastic::Solve(CRegion &Region,double DeltaT)
{
	//UpdateDs(Region);
	//UpdateS(Region,DeltaT);

	//CSPHPt * PtiPtr,* PtjPtr;
	//CKnl * KnlPtr;

	//double normx;
	//double normy;
	//double normz;


	//for (size_t i=0;i<Region._PtPairList.size();i++)
	//{
	//	if(Region._PtPairList[i]._Type==enSPHPtPair)
	//	{
	//		PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
	//		PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
	//		KnlPtr=&Region._KnlList[i];

	//		normx=((PtiPtr->_sigmaxx+PtjPtr->_sigmaxx)*KnlPtr->_Wx+(PtiPtr->_sigmaxy+PtjPtr->_sigmaxy)*KnlPtr->_Wy+(PtiPtr->_sigmaxz+PtjPtr->_sigmaxz)*KnlPtr->_Wz)/(PtiPtr->_rho*PtjPtr->_rho);
	//		normy=((PtiPtr->_sigmaxy+PtjPtr->_sigmaxy)*KnlPtr->_Wx+(PtiPtr->_sigmayy+PtjPtr->_sigmayy)*KnlPtr->_Wy+(PtiPtr->_sigmayz+PtjPtr->_sigmayz)*KnlPtr->_Wz)/(PtiPtr->_rho*PtjPtr->_rho);
	//		normz=((PtiPtr->_sigmaxz+PtjPtr->_sigmaxz)*KnlPtr->_Wx+(PtiPtr->_sigmayz+PtjPtr->_sigmayz)*KnlPtr->_Wy+(PtiPtr->_sigmazz+PtjPtr->_sigmazz)*KnlPtr->_Wz)/(PtiPtr->_rho*PtjPtr->_rho);

	//		PtiPtr->_du+=PtjPtr->_m*normx;
	//		PtiPtr->_dv+=PtjPtr->_m*normy;
	//		PtiPtr->_dw+=PtjPtr->_m*normz;
	//		//PtiPtr->_de+=???
	//		if(PtiPtr!=PtjPtr)
	//		{
	//			PtjPtr->_du-=PtiPtr->_m*normx;
	//			PtjPtr->_dv-=PtiPtr->_m*normy;
	//			PtjPtr->_dw-=PtiPtr->_m*normz;
	//			//PtiPtr->_de+=???
	//		}
	//	}
	//}
}
