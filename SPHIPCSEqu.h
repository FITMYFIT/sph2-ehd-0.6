/******************************************************************************
!
! Purpose:to solve the SPHIPCS equations  
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHIPCSEQU,v 2.0 (date)03.09 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHIPCSEQU_H_
#define _SPHIPCSEQU_H_

#include "SPHPt.h"
#include "Region.h"
#include "SPHIPCS.h"

class CSPHIPCSEqu
{
public:
	CSPHIPCSEqu();

	~CSPHIPCSEqu();

	void Solve(CRegion & Region,double DeltaT);

private:

	void Init(CRegion & Region);
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
