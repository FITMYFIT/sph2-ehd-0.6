/******************************************************************************
!
! Purpose:to solve the equation of interaction force  
!
! Description:the basic part of N-S equation
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHINTFORCE,v 2.0 (date)03.11 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _CONTINUITYEQU_H_
#define _CONTINUITYEQU_H_

#include "Region.h"
#include "SPHPt.h"
#include "BasePt.h"
#include "CSPM.h"
#include "GetKnlList.h"

class CContinuityEqu
{
public:

	CContinuityEqu();

	~CContinuityEqu();

	void Solve(CRegion &Region,unsigned int TimeSteps);

private:
	CCSPM _CSPM;

	CGetKnlList _GetKnlListC;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
