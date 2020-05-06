/******************************************************************************
!
! Purpose:to solve some special problems
!
! Description:the same as SectionSPH in content
!
! Notes:
!
!******************************************************************************
!
!$Id: SECTIONNULL,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SECTIONNULL_H_
#define _SECTIONNULL_H_

#include "Section.h"

class CSectionNULL:public CSection
{
public:
	
	CSectionNULL(unsigned int SECID=0);
	
	~CSectionNULL();
	
	double _CSLH;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
