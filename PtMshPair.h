/******************************************************************************
!
! Purpose:
!
! Description:the basic type of the PtPairList
!
! Notes:
!
!******************************************************************************
!
!$Id: PTPAIR,v 2.0 (date)03.08 CHEN Fuzhen
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _PTMSHPAIR_H_
#define _PTMSHPAIR_H_

#include "BasePt.h"
#include "Mesh.h"


class CPtMshPair
{
public:
	CPtMshPair();

	~CPtMshPair();

	CBasePt * _PtiPtr;

	CMesh * _MshjPtr;

	double _driac2;

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
