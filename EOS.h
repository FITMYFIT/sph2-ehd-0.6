/******************************************************************************
!
! Purpose:to define the basic type of the equation of state
!
! Description:the equation of state
!
! Notes:
!
!******************************************************************************
!
!$Id: EOS,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _EOS_H_
#define _EOS_H_

enum enEOSTYPE {enIdealGas,enWeaklyCompress,enJWL,enIPCS};

class CEOS
{
public:
	CEOS(enEOSTYPE Type,unsigned int EOSID=0);

	~CEOS();

	enEOSTYPE Type();

	unsigned int _EOSID;

private:
	enEOSTYPE _Type;
};

#endif 

/******************************************************************************

!******************************************************************************/
