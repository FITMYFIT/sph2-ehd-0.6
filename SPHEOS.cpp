#include "SPHEOS.h"

CSPHEOS::CSPHEOS()
{
}

CSPHEOS::~CSPHEOS()
{
}

void CSPHEOS::Solve(CRegion &Region)
{
	unsigned int jj;
	CPart * PartPtr;
	CSPHPt * SPHPtPtr;
	CBasePt * BasePtPtr;
	CEOS * EOSPtr;
	unsigned int ID,PID;

	for(size_t i=0;i!=Region._PtList.size();++i)
	{
		BasePtPtr=&Region._PtList[i];

		//if(BasePtPtr->_Type==enSPHPt||BasePtPtr->_Type==enNULLPt)
		{
			PID=BasePtPtr->_PID;

			EOSPtr=Region._EOSList[Region._PartList[PID-1]._EOSID-1];
			
			if(EOSPtr->Type()==enIdealGas)
			{
				((CIdealGas *)EOSPtr)->EOS(BasePtPtr->_rho,BasePtPtr->_e,BasePtPtr->_p,BasePtPtr->_Cs);
			}
			else if(EOSPtr->Type()==enWeaklyCompress)
			{
				((CWeaklyCompress *)EOSPtr)->EOS(BasePtPtr->_rho0,BasePtPtr->_rho,BasePtPtr->_p);
				//((CWeaklyCompress *)EOSPtr)->EOS2(Region._ControlSPH._Cs,BasePtPtr->_rho,BasePtPtr->_p);//p=cs*cs*p
			}
		}
	}
}
