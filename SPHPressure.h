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

#ifndef _SPHPRESSURE_H_
#define _SPHPRESSURE_H_

#include "Region.h"
#include "SPHPt.h"
#include "GetKnlList.h"

class CSPHPressure
{
public:
	CSPHPressure();

	~CSPHPressure();

	void Solve(CRegion &Region);

private:
  
	CGetKnlList _GetKnlListP;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
