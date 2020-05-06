/******************************************************************************
!
! Purpose:to solve the equation of state
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHEOS,v 2.0 (date)03.09 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHEOS_H_
#define _SPHEOS_H_

#include "IdealGas.h"
#include "WeaklyCompress.h"
#include "Region.h"
#include "SPHPt.h"
#include "Part.h"
#include "EOS.h"
#include <vector>

class CSPHEOS
{
public:
	CSPHEOS();

	~CSPHEOS();

	void Solve(CRegion &Region);

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
