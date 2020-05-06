/******************************************************************************
!
! Purpose:consindering the shearingforce  
!
! Description:Jaumann ratio
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHSHEARINGFORCE,v 2.0 (date)03.17 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHSHEARINGFORCE_H_
#define _SPHSHEARINGFORCE_H_

#include "SPHPt.h"
#include "Region.h"
#include <cmath>

class CSPHShearingForce
{
public:
	CSPHShearingForce();

	~CSPHShearingForce();

	void Solve(CRegion &Region,double DeltaT);

private:
	void UpdateDs(CRegion &Region);

	void UpdateS(CRegion &Region,double DeltaT);
};


#endif

/******************************************************************************
!
!
!******************************************************************************/
