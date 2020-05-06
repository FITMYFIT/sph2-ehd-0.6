#include "SPHIPCS.h"

CSPHIPCS::CSPHIPCS(unsigned int EOSID,double c,double somega,double bomega,unsigned int SurfaceBool,double SurfaceCoeff,unsigned int ShearingBool,double ViscocityEtas)
:CEOS(enIPCS,EOSID),_c(c),_somega(somega),_bomega(bomega),_SurfaceBool(SurfaceBool),_SurfaceCoeff(SurfaceCoeff),_ShearingBool(ShearingBool),
_ViscocityEtas(ViscocityEtas)
{
}

CSPHIPCS::~CSPHIPCS()
{
}

void CSPHIPCS::GetDeltaP(CSPHPt *PtiPtr,CSPHPt *PtjPtr,CKnl * KnlPtr,double DeltaT)
{
	double uij;
	double vij;
	double wij;
	double temp;

	uij=PtiPtr->_u-PtjPtr->_u;
	vij=PtiPtr->_v-PtjPtr->_v;
	wij=PtiPtr->_w-PtjPtr->_w;
	temp=uij*KnlPtr->_Wx+vij*KnlPtr->_Wy+wij*KnlPtr->_Wz;

	PtiPtr->_deltap1+=_c*_c*DeltaT*PtjPtr->_m*temp;
	if(PtiPtr!=PtjPtr)
	{
		PtjPtr->_deltap1+=_c*_c*DeltaT*PtiPtr->_m*temp;
	}
}

void CSPHIPCS::ReInit(CSPHPt * SPHPtPtr)
{
	SPHPtPtr->_deltap=SPHPtPtr->_deltap1;

	SPHPtPtr->_deltap1=0.0;
	SPHPtPtr->_deltau=0.0;
	SPHPtPtr->_deltav=0.0;
	SPHPtPtr->_deltaw=0.0;
}

void CSPHIPCS::GetDeltaV(CSPHPt *PtiPtr, CSPHPt *PtjPtr, CKnl *KnlPtr, double DeltaT)
{
	PtiPtr->_deltau-=(DeltaT/PtiPtr->_rho)*(PtjPtr->_m/PtjPtr->_rho)*(PtjPtr->_deltap-PtiPtr->_deltap)*KnlPtr->_Wx;
	PtiPtr->_deltav-=(DeltaT/PtiPtr->_rho)*(PtjPtr->_m/PtjPtr->_rho)*(PtjPtr->_deltap-PtiPtr->_deltap)*KnlPtr->_Wy;
	PtiPtr->_deltaw-=(DeltaT/PtiPtr->_rho)*(PtjPtr->_m/PtjPtr->_rho)*(PtjPtr->_deltap-PtiPtr->_deltap)*KnlPtr->_Wz;

	if(PtiPtr!=PtjPtr)
	{
		PtjPtr->_deltau-=(DeltaT/PtjPtr->_rho)*(PtiPtr->_m/PtiPtr->_rho)*(PtjPtr->_deltap-PtiPtr->_deltap)*KnlPtr->_Wx;
		PtjPtr->_deltav-=(DeltaT/PtjPtr->_rho)*(PtiPtr->_m/PtiPtr->_rho)*(PtjPtr->_deltap-PtiPtr->_deltap)*KnlPtr->_Wy;
		PtjPtr->_deltaw-=(DeltaT/PtjPtr->_rho)*(PtiPtr->_m/PtiPtr->_rho)*(PtjPtr->_deltap-PtiPtr->_deltap)*KnlPtr->_Wz;
	}
}

void CSPHIPCS::UpdateV(CSPHPt *SPHPtPtr)
{
	SPHPtPtr->_u+=_bomega*SPHPtPtr->_deltau;
	SPHPtPtr->_v+=_bomega*SPHPtPtr->_deltav;
	SPHPtPtr->_w+=_bomega*SPHPtPtr->_deltaw;
}

void CSPHIPCS::GetSMax(CSPHPt *SPHPtPtr,double &SMax)
{
	double S;

	if(SPHPtPtr->_deltap1>=SPHPtPtr->_deltap)
	{
		S=SPHPtPtr->_deltap1-SPHPtPtr->_deltap;
	}
	else
	{
		S=SPHPtPtr->_deltap-SPHPtPtr->_deltap1;
	}
	if(S>SMax)
	{
		SMax=S;
	}
}

void CSPHIPCS::UpdateP(CSPHPt *SPHPtPtr)
{
	SPHPtPtr->_p+=_somega*SPHPtPtr->_deltap;
}
