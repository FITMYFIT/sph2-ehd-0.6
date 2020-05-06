#include "DeltaT.h"

CDeltaT::CDeltaT()
:_DeltaT(0.0),_DeltaT1(0.0)
{
}

CDeltaT::~CDeltaT()
{
}

void CDeltaT::GetDeltaT(CRegion & Region)
{
	_DeltaT=_DeltaT1=Region._ControlSPH._DeltaT;

    //unsigned int jj;
	//CSPHPt * SPHPtPtr;
	//CPart * PartPtr;
	//double DeltaTT;//the middle parametric of DeltaT

	//_DeltaT=_DeltaT1;
	//jj=Region._CalList[0];
 //   SPHPtPtr=((CSPHPt *)Region._PtList[Region._CalList[0]]);
	//_DeltaT1=Region._ControlSPH._DeltaTCoeff*SPHPtPtr->_h/SPHPtPtr->_c;
	//
	//for (size_t i=0;i<Region._PartList.size();i++)
	//{
	//	PartPtr=(CPart *)Region._PartList[i];

	//	if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
	//		||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
	//	{
	//		for(size_t j=0;j<PartPtr->_PartCalList.size();j++)
	//		{
	//			jj=PartPtr->_PartCalList[j];
	//			SPHPtPtr=(CSPHPt *)PartPtr->_PartPtList[jj];
	//			DeltaTT=Region._ControlSPH._DeltaTCoeff*SPHPtPtr->_h/SPHPtPtr->_c;

	//			if(DeltaTT<_DeltaT1)
	//				_DeltaT1=DeltaTT;

	//			if(_DeltaT1<=0.0000000001)
	//				_DeltaT1=0.0000000001;
	//			if(_DeltaT1>=0.00000001)
	//				_DeltaT1=0.00000001;
	//		}
	//	}
	//}
}
