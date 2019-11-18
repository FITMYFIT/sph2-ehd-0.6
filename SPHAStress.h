/******************************************************************************
!
! Purpose:to remove the tensor instability
!
! Description:artifical stress derived by Monaghan(2000) and extened by Gray(2001)
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHASTRESS,v 2.0 (date)03.11 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHASTRESS_H_
#define _SPHASTRESS_H_

#include "Region.h"
#include "SPHPt.h"
#include "GetKnlList.h"

class CSPHAStress
{
public:
	CSPHAStress();

	~CSPHAStress();

	void Solve(CRegion &Region);

private:
	CGetKnlList _GetKnlList;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
