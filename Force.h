/******************************************************************************
!
! Purpose:to define the basic type of force
!
! Description:the basic type of the _ExtForceList
!
! Notes:
!
!******************************************************************************
!
!$Id: FORCE,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _FORCE_H_
#define _FORCE_H_
#include<vector>

class CForce
{
public:
	CForce(size_t PID=0);
	
	~CForce();

	size_t _PID;
	
	double _fx;
	
	double _fy;

	double _fz;
private:	
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
