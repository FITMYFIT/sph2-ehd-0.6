/******************************************************************************
!
! Purpose:to define a type of equation of state
!
! Description:the equation of the state of IdealGas
!
! Notes:
!
!******************************************************************************
!
!$Id: IDEALGAS,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _IDEALGAS_H_
#define _IDEALGAS_H_


#include "EOS.h"
#include <cmath>
class CIdealGas:public CEOS
{
public:
	CIdealGas(double Cp,double Mw,double R,unsigned int EOSID=0);

	~CIdealGas();	
	
	double _Cp;

	double _Mw;//mass of molecule

	double _R;

    //the equation of state related to energy
	void EOS(double rho,double ee,double &p,double &c);

  //the equation of state related to temperature
	void EOS2(double rho,double T,double &p,double &c);

private:
	double _gama;

	double _gama1;//_gama-1.0

	double _Rg;//_R/_Mw

	double _Cv;//_Cp-_Rg

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
