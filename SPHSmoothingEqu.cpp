#include "SPHSmoothingEqu.h"

CSPHSmoothingEqu::CSPHSmoothingEqu()
{
}

CSPHSmoothingEqu::~CSPHSmoothingEqu()
{
}

void CSPHSmoothingEqu::Solve(CRegion &Region)
{
	unsigned int jj;
	CSPHPt * SPHPtPtr;
	CPart * PartPtr;

	for (size_t i=0;i<Region._PartList.size();i++)
	{
		PartPtr=&Region._PartList[i];

		if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
			||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
		{
			for(size_t j=0;j<PartPtr->_PartCalList.size();j++)
			{
				jj=PartPtr->_PartCalList[j];
				SPHPtPtr=PartPtr->_PartPtList[jj];
				SPHPtPtr->_dh=-0.5*SPHPtPtr->_h*SPHPtPtr->_drho/SPHPtPtr->_rho;
			}
		}
	}
}
