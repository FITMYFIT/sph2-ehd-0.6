/******************************************************************************
!
! Purpose: to define a type of particles
!
! Description:the type of SPH particles
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHPT,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHINIT_H_
#define _SPHINIT_H_

#include "BasePt.h"
#include "Region.h"
#include <vector>

class CSPHInit
{
public:

	CSPHInit();

	~CSPHInit();

	void Solve(CRegion &Region,unsigned int TimeSteps);//≥ı ºªØ

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
