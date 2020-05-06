/******************************************************************************
!
! Purpose:to define the node parameters
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: CBaseGridNode,v 2.0 (date)03.11. CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _BASEGRIDNODE_H_
#define _BASEGRIDNODE_H_

class CBaseGridNode
{
public:
	
	CBaseGridNode(unsigned int ID=0);
	
	~CBaseGridNode();

	unsigned int _ID;
	
	double _x;
	
	double _y;

	double _z;

	double _C;

	double _Cs;

	double _div;

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
