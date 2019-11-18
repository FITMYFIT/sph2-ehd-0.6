/******************************************************************************
!
! Purpose:to get the timestep size 
!
! Description:according to  CFL condition
!
! Notes:
!
!******************************************************************************
!
!$Id: DELTAT,v 2.0 (date)03.09 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _DELTAT_H_
#define _DELTAT_H_

#include "Region.h"
#include "SPHPt.h"

class CDeltaT
{
public:
	CDeltaT();

	~CDeltaT();

	double _DeltaT;//_DeltaT(n)

	double _DeltaT1;//_DeltaT(n+1)

	void GetDeltaT(CRegion & Region);

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
