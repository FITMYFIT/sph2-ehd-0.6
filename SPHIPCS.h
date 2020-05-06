/******************************************************************************
!
! Purpose:to solve the equation of improved pressure correction scheme(H.S.Fang2009)
!
! Description:combination of the equation of state method and the global pressure Poisson method
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHIPCS,v 2.0 (date)03.09 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHIPCS_H_
#define _SPHIPCS_H_

#include "Region.h"
#include "SPHPt.h"
#include "EOS.h"

class CSPHIPCS:public CEOS
{
public:
	CSPHIPCS(unsigned int EOSID,double c,double somega,double bomega,unsigned int SurfaceBool,double SurfaceCoeff,unsigned int ShearingBool,double ViscocityEtas);

	~CSPHIPCS();

	unsigned int _SurfaceBool;

	double _SurfaceCoeff;

	unsigned int _ShearingBool;

	double _ViscocityEtas;

	void GetDeltaP(CSPHPt *PtiPtr,CSPHPt *PtjPtr,CKnl * KnlPtr,double DeltaT);

	void ReInit(CSPHPt * SPHPtPtr);

	void GetDeltaV(CSPHPt *PtiPtr,CSPHPt *PtjPtr,CKnl * KnlPtr,double DeltaT);

	void UpdateV(CSPHPt * SPHPtPtr);

	void GetSMax(CSPHPt * SPHPtPtr,double &SMax);
	
	void UpdateP(CSPHPt * SPHPtPtr);

private:
	double _c;

	double _somega;

	double _bomega;

	double _DeltaT;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
