/******************************************************************************
!
! Purpose:to transform kinetic energy to thermal energy and avoid non-physical penetration
!
! Description:artifical viscous force derived by Monaghan(1985)
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHAVFORCE,v 2.0 (date)03.11 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHAVFORCE_H_
#define _SPHAVFORCE_H_

#include "Region.h"
#include "SPHPt.h"
#include <cmath>

class CSPHAVForce
{
public:
	CSPHAVForce();

	~CSPHAVForce();

	void Solve(CRegion &Region);
	void Solve2(CRegion &Region);//Monaghan1997型人工粘性
private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
