#include "WeaklyCompress.h"

CWeaklyCompress::CWeaklyCompress(double P0,double gamma,double ExtP,unsigned int EOSID)
:CEOS(enWeaklyCompress,EOSID),_P0(P0),_gamma(gamma),_ExtP(ExtP)
{
}

CWeaklyCompress::~CWeaklyCompress()
{
}

void CWeaklyCompress::EOS(double rho0,double rho,double &p)
{
	double rhol;//rho/rho0
	rhol=rho/rho0;
	p=_P0*(pow(rhol,_gamma)-1)+_ExtP;
}

void CWeaklyCompress::EOS2( double cs,double rho,double &p )
{

	p=cs*cs*rho;

}
