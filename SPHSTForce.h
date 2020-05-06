/******************************************************************************
!
! Purpose:to solve the problems which consider the surface tension model
!
! Description:Van der Waals cohensive force applied by Nugent and Posch(2000)
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHSURFACEFORCE4,v 2.0 (date)03.11 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHSTFORCE_H_
#define _SPHSTFORCE_H_

#include "Region.h"
#include "BasePt.h"
#include "SPHPt.h"
#include "CSPM.h"
#include "Operator.h"
#include <vector>

class CSPHSTForce
{
public:
	CSPHSTForce();

	~CSPHSTForce();

	void Solve(CRegion &Region);

private:

	CCSPM _CSPMST;

	COperator _OperatorST;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
