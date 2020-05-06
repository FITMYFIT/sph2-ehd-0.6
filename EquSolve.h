/******************************************************************************
!
! Purpose:迭代求解方程
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: EQUSOLVE,v 2.0 (date)03.18 Richard LIU $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _EQUSOLVE_H_
#define _EQUSOLVE_H_

#include "SPHPt.h"
#include "Region.h"
#include "Operator.h"
#include "Matrix.h"
#include <cmath>
#include <vector>
#include <iostream>

class CEquSolve
{
public:

	CEquSolve();

	~CEquSolve();

	void BICGSolve(const std::vector<CMatrix> &A,std::vector<double> &x,const std::vector<double> &b,unsigned int MaxIter,double RTC);

private:

  COperator _Operator;

	unsigned int RTCCtrl(double rNorm,double bNorm,double RTC);//精度控制，达到精度后退出


};

#endif

/******************************************************************************
!
!
!******************************************************************************/
