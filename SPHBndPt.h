/******************************************************************************
!
! Purpose:to define a type of particle
!
! Description:the type of SPH boundary particles
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHBNDPT,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHBNDPT_H_
#define _SPHBNDPT_H_

#include "BasePt.h"

class CSPHBndPt : public CBasePt
{
public:
	CSPHBndPt(unsigned int ID,unsigned int PartID);
	
	~CSPHBndPt();
	

	double _nx;//the direction

	double _ny;

	double _nz;

	double _area;
private:	
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
