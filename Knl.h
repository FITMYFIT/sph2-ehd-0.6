/******************************************************************************
!
! Purpose:to define the basic type of kernel
!
! Description:contain some parameters of kernel
!
! Notes:
!
!******************************************************************************
!
!$Id: KNL,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _KNL_H_
#define _KNL_H_

class CKnl
{
public:
	CKnl();
	
	~CKnl();

	double _W;

	double _Wx;

	double _Wy;

	double _Wz;

	double _W2;

	double _Ww;

	double _Wj;//用于变光滑长度

	double _Wxj;

	double _Wyj;

	double _Wzj;

	double _W2j;

	double _Wwj;

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
