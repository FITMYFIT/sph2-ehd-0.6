/******************************************************************************
!
! Purpose:to control some command related to SPH particles
!
! Description:a type of section
!
! Notes:
!
!******************************************************************************
!
!$Id: SECTIONSPH,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SECTIONSPH_H_
#define _SECTIONSPH_H_

#include "Section.h"

class CSectionSPH:public CSection
{
public:
	
	CSectionSPH(unsigned int SECID=0);
	
	~CSectionSPH();
	
	double _CSLH;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
