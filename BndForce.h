/******************************************************************************
!
! Purpose:to solve the problems considering the boundary force
!
! Description:the boundary force
!
! Notes:
!
!******************************************************************************
!
!$Id: BNDFORCE,v 2.0 (date)03.14 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _BNDFORCE_H_
#define _BNDFORCE_H_

#include "Region.h"
#include "SPHPt.h"
#include "Operator.h"

class CBndForce
{
public:
	CBndForce();

	~CBndForce();

	void Solve(CRegion &Region);

private:
	COperator _OperatorB;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
