/******************************************************************************
!
! Purpose:to define the basic type of the section
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: SECTION,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SECTION_H_
#define _SECTION_H_

enum enSECTIONTYPE {enSectionSPH,enSectionNULL};

class CSection
{
public:
	CSection(enSECTIONTYPE Type,unsigned int SECID=0);
	
	~CSection();

	enSECTIONTYPE Type();

	unsigned int _SECID;

private:
	enSECTIONTYPE _Type;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
