#include "IdealGas.h"
#include <vector>

CIdealGas::CIdealGas(double Cp,double Mw,double R,unsigned int EOSID)
:CEOS(enIdealGas,EOSID),_Cp(Cp),_Mw(Mw),_R(R)
{
    _Rg=_R/_Mw;
	_Cv=_Cp-_Rg;
	_gama1=_Rg/_Cv;
    _gama=_gama1+1.0;
}

CIdealGas::~CIdealGas()
{
}

void CIdealGas::EOS(double rho,double ee,double &p,double &c)
{
	p=_gama1*rho*ee;
	c=sqrt(std::max(_gama*_gama1*ee,0.0));
}

void CIdealGas::EOS2(double rho,double T,double &p,double &c)
{
	p=_gama1*rho*_Cv*T;
	c=sqrt(std::max(_gama1*(_gama1-1)*_Cv*T,0.0));
}
