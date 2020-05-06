/******************************************************************************
!
! Purpose:consindering the shearingforce  
!
! Description:viscoelastic model Maxwell model(M.Ellero2002¡¢2005 and A.Rafiee2007)
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHVISCOELASTIC,v 2.0 (date)03.18 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHVISCOELASTIC_H_
#define _SPHVISCOELASTIC_H_

#include "SPHPt.h"
#include "Region.h"
#include <cmath>

class CSPHViscoelastic
{
public:
	CSPHViscoelastic();

	~CSPHViscoelastic();

	void Solve(CRegion &Region,double DeltaT);

	void UpdateDs(CRegion &Region);

	void UpdateS(CRegion &Region,double DeltaT);

private:

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
