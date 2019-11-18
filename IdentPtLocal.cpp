#include "IdentPtLocal.h"

CIdentPtLocal::CIdentPtLocal()
{

}

CIdentPtLocal::~CIdentPtLocal()
{

}


void CIdentPtLocal::IdentPtLocal(CRegion & Region)
{
	if(Region._ControlSPH._CSPMIflag2==0)
	{
		_CSPMI.GetCSPMGradCorctCoef(Region);
	}

	CBasePt * BasePtPtr;

	unsigned int ID;

	double Axx,Axy,Ayx,Ayy;

	if(Region._ControlSPH._CSPMIflag2==0)
	{
		_CSPMI.GetCSPMGradCorctCoef(Region);
	}


	for(unsigned int i=0;i!=Region._CalList.size();++i)
	{
		ID=Region._CalList[i];

		BasePtPtr=&Region._PtList[ID];

		Axx=-BasePtPtr->_CSPMAxx;
		Ayx=-BasePtPtr->_CSPMAyx;
		Axy=-BasePtPtr->_CSPMAxy;
		Ayy=-BasePtPtr->_CSPMAyy;

		//_OperaterI.Reve2ndMat(Axx,Axy,Ayx,Ayy,&Axx,&Axy,&Ayx,&Ayy);//Çó³öB¾ØÕó£¬S Marrone 2010 Equ1

		BasePtPtr->_CSPMminLambda=_OperaterI.MinEig2ndMat(Axx,Axy,Ayx,Ayy);
	}

}