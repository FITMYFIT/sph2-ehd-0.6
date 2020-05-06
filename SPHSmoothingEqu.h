/******************************************************************************
!
! Purpose:to solve the equation of changing smoothing length 
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHSMOOTHINGEQU,v 2.0 (date)03.09 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHSMOOTHINGEQU_H_
#define _SPHSMOOTHINGEQU_H_

#include "Region.h"
#include "SPHPt.h"

class CSPHSmoothingEqu
{
public:
	CSPHSmoothingEqu();

	~CSPHSmoothingEqu();

	void Solve(CRegion &Region);

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
