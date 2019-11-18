/******************************************************************************
!
! Purpose:to get the kernel list
!
! Description:the kernel type is B3SKn
!
! Notes:
!
!******************************************************************************
!
!$Id: GETKNLLIST,v 2.0 (date)03.09 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _GETKNLLIST_H_
#define _GETKNLLIST_H_

#include "Region.h"
#include "BasePt.h"
#include "Knl.h"
#include "SPHPt.h"
#include "CSPM.h"
#include "Operator.h"
#include "Mesh.h"
#include <cmath>

enum enKNLTYPE {enB3SKnl, enB4NSKnl,enB4SPKnl,enB5SPKnl,enBGKnl};


class CGetKnlList
{
public:

	CGetKnlList(enKNLTYPE Type=enB3SKnl);

	~CGetKnlList();

	enKNLTYPE _Type;

	void GetKnlList(CRegion &Region);
	void GetMshKnlList(CRegion &Region);

	void GetCrctKnlGrad(CRegion & Region,CBasePt * PtPtr,CKnl * KnlPtr,double & wcx,double & wcy);//ÐÞÕýºËº¯ÊýÌÝ¶È

	void ClearKnlList(CRegion &Region);

	void GetW0(double distance0,double h,double & W);

private:
  double _norm;

	CCSPM _CSPMG;

  COperator _OperatorG;

	void GetW(double x,double y,double h,double distance2,double & W,double & Wx,double & Wy,double & Ww);

	void GetW2(double x,double y,double h,double distance2,double & W,double & Wx,double & Wy,double & W2);
};
#endif

/******************************************************************************
!
!
!******************************************************************************/
