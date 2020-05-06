/******************************************************************************
!
! Purpose:to update the parameters of particles
!
! Description:leap-frog method
!
! Notes:
!
!******************************************************************************
!
!$Id: LEAPFROG,v 2.0 (date)03.11 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _UPDATEPOSITION_H_
#define _UPDATEPOSITION_H_

#include "Region.h"
#include "DeltaT.h"
#include "SPHPt.h"
#include "BasePt.h"
#include "Knl.h"
#include <cmath>

class CUpdatePosition
{
public:
	CUpdatePosition();

	~CUpdatePosition();

	void LeapFrogUpdate(CRegion & Region, double DeltaT,double DeltaT1,double TimeSteps);

	void ImplicitUpdate(CRegion & Region);

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
