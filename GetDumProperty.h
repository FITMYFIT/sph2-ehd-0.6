/******************************************************************************
!
! Purpose:X Y Hu2012dummy particle method to implement boundary condition
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

#ifndef _GETDUMPROPERTY_H_
#define _GETDUMPROPERTY_H_

#include "Region.h"
#include "SPHPt.h"
#include "Operator.h"

#include <cmath>

class CGetDumProperty
{
public:
	CGetDumProperty();

	~CGetDumProperty();

	void Solve(CRegion &Region);

private:
	COperator _OperatorGD;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
