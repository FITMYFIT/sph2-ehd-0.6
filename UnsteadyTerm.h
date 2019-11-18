/******************************************************************************
!
! Purpose:consindering the shearingforce
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHVISCOCITY,v 2.0 (date)03.18 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _UNSTEADYTERM_H_
#define _UNSTEADYTERM_H_

#include "SPHPt.h"
#include "Region.h"
#include "Matrix.h"
#include <cmath>

class CUnsteadyTerm
{
public:

	CUnsteadyTerm();

	~CUnsteadyTerm();

	void Solve(CRegion &Region);

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
