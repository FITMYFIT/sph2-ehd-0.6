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

#ifndef _PTPAIR_H_
#define _PTPAIR_H_

#include "BasePt.h"

enum enPTPAIRTYPE {enSPHPtPair,enSPHNULLPtPair, enSPHBndPtPair,enSPHGhostPtPair,
                   enSPHDumPtPair,enSPHEHDDumPtPair,enSPHEHDBndPtPair,enEHDBndDumPtPair};

class CPtPair
{
public:
	CPtPair();

	CPtPair(enPTPAIRTYPE Type, CBasePt * PtiPtr, CBasePt * PtjPtr,double driac2);
	
	~CPtPair();
	
	enPTPAIRTYPE _Type;
	
	CBasePt * _PtiPtr;
	
	CBasePt * _PtjPtr;

	double _driac2;

  double _xij,_yij,_zij;//for simpicity when periodic boundary condition is implemented

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
