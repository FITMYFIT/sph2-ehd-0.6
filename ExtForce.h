/******************************************************************************
!
! Purpose:
!
! Description:the extral force
!
! Notes:
!
!******************************************************************************
!
!$Id: EXTFORCE,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _EXTFORCE_H_
#define _EXTFORCE_H_

#include "Region.h"
#include "SPHPt.h"

#include <cmath>

class CExtForce
{
public:
	CExtForce();

	~CExtForce();

	void Solve(CRegion &Region);
	void Solve2(CRegion &Region,double TimeSteps);//加入damp的外力施加，XY Hu2012 Equ13

private:

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
