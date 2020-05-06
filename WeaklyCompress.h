/******************************************************************************
!
! Purpose:define a type of equation of state
!
! Description:weakly compressibility equation applied by Monaghan(1994)
!
! Notes:
!
!******************************************************************************
!
!$Id: WEAKLYCOMPRESS,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _WEAKLYCOMPRESS_H_
#define _WEAKLYCOMPRESS_H_

#include "EOS.h"
#include <cmath>

class CWeaklyCompress:public CEOS
{
public:
	CWeaklyCompress(double P0,double gamma,double ExtP,unsigned int EOSID);

	~CWeaklyCompress();

	unsigned int _SurfaceBool;
	unsigned int _SurfaceCoeff;

	void EOS(double rho0,double rho,double &p);
	void EOS2(double cs,double rho,double &p);//morris, P=cs^2*rho

private:
	double _P0;

	double _gamma;

	double _ExtP;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
