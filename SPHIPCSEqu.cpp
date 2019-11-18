#include "SPHIPCSEqu.h"

CSPHIPCSEqu::CSPHIPCSEqu()
{
}

CSPHIPCSEqu::~CSPHIPCSEqu()
{
}

void CSPHIPCSEqu::Init(CRegion &Region)
{
	CSPHPt * SPHPtPtr;
	CPart * PartPtr;

	for (size_t i=0;i<Region._PartList.size();i++)
	{
		PartPtr=&Region._PartList[i];

		if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
			||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
		{
			for(size_t j=0;j<PartPtr->_PartPtList.size();j++)
			{
				SPHPtPtr=PartPtr->_PartPtList[j];

				SPHPtPtr->_deltap=0.0;
				SPHPtPtr->_deltap1=0.0;
				SPHPtPtr->_deltau=0.0;
				SPHPtPtr->_deltav=0.0;
				SPHPtPtr->_deltaw=0.0;
			}
		}
	}
}

void CSPHIPCSEqu::Solve(CRegion &Region,double DeltaT)
{
	CPart * PartPtr;
	CSPHPt * SPHPtPtr;
	CEOS * EOSPtr;
	CSPHPt *PtiPtr,*PtjPtr;
	CKnl *KnlPtr;
	CEOS *EOSiPtr,*EOSjPtr;

	double SMax;

	Init(Region);

	for (size_t i=0;i<Region._PtPairList.size();i++)
	{
		if(Region._PtPairList[i]._Type==enSPHPtPair)
		{
			PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
			PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
			KnlPtr=&Region._KnlList[i];

			EOSiPtr=Region._EOSList[Region._PartList[PtiPtr->_PID-1]._EOSID-1];
			EOSjPtr=Region._EOSList[Region._PartList[PtjPtr->_PID-1]._EOSID-1];
			if(EOSiPtr==EOSjPtr)
			{
				if(EOSiPtr->Type()==enIPCS)
				{
					((CSPHIPCS *)EOSiPtr)->GetDeltaP(PtiPtr,PtjPtr,KnlPtr,DeltaT);
				}
			}
		}
	}

	do
	{
		SMax=0.0000;
		for (size_t i=0;i<Region._PartList.size();i++)
		{
			PartPtr=&Region._PartList[i];
			if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
				||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
			{
				for(size_t j=0;j<PartPtr->_PartPtList.size();j++)
				{
					SPHPtPtr=PartPtr->_PartPtList[j];
					EOSPtr=Region._EOSList[PartPtr->_EOSID-1];
					
					if(EOSPtr->Type()==enIPCS)
					{
						((CSPHIPCS *)EOSPtr)->ReInit(SPHPtPtr);
					}
				}
			}
		}

		for (size_t i=0;i<Region._PtPairList.size();i++)
		{
			if(Region._PtPairList[i]._Type==enSPHPtPair)
			{
				PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
				PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
				KnlPtr=&Region._KnlList[i];

				EOSiPtr=Region._EOSList[Region._PartList[PtiPtr->_PID-1]._EOSID-1];
				EOSjPtr=Region._EOSList[Region._PartList[PtjPtr->_PID-1]._EOSID-1];
				if(EOSiPtr==EOSjPtr)
				{
					if(EOSiPtr->Type()==enIPCS)
					{
						((CSPHIPCS *)EOSiPtr)->GetDeltaV(PtiPtr,PtjPtr,KnlPtr,DeltaT);
					}
				}
			}
		}

		for (size_t i=0;i<Region._PartList.size();i++)
		{
			PartPtr=&Region._PartList[i];
			if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
				||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
			{
				for(size_t j=0;j<PartPtr->_PartPtList.size();j++)
				{
					SPHPtPtr=PartPtr->_PartPtList[j];
					EOSPtr=Region._EOSList[PartPtr->_EOSID-1];
					
					if(EOSPtr->Type()==enIPCS)
					{
						((CSPHIPCS *)EOSPtr)->UpdateV(SPHPtPtr);
					}
				}
			}
		}

		for (size_t i=0;i<Region._PtPairList.size();i++)
		{
			if(Region._PtPairList[i]._Type==enSPHPtPair)
			{
				PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
				PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
				KnlPtr=&Region._KnlList[i];

				EOSiPtr=Region._EOSList[Region._PartList[PtiPtr->_PID-1]._EOSID-1];
				EOSjPtr=Region._EOSList[Region._PartList[PtjPtr->_PID-1]._EOSID-1];
				if(EOSiPtr==EOSjPtr)
				{
					if(EOSiPtr->Type()==enIPCS)
					{
						((CSPHIPCS *)EOSiPtr)->GetDeltaP(PtiPtr,PtjPtr,KnlPtr,DeltaT);
					}
				}
			}
		}

		for (size_t i=0;i<Region._PartList.size();i++)
		{
			PartPtr=&Region._PartList[i];
			if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
				||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
			{
				for(size_t j=0;j<PartPtr->_PartPtList.size();j++)
				{
					SPHPtPtr=PartPtr->_PartPtList[j];
					EOSPtr=Region._EOSList[PartPtr->_EOSID-1];
					
					if(EOSPtr->Type()==enIPCS)
					{
						((CSPHIPCS *)EOSPtr)->GetSMax(SPHPtPtr,SMax);
					}
				}
			}
		}
	}while(SMax>100.0);

	for (size_t i=0;i<Region._PartList.size();i++)
	{
		PartPtr=&Region._PartList[i];
		if(Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionSPH
			||Region._SectionList[PartPtr->_SECID-1]->Type()==enSectionNULL)
		{
			for(size_t j=0;j<PartPtr->_PartPtList.size();j++)
			{
				SPHPtPtr=PartPtr->_PartPtList[j];
				EOSPtr=Region._EOSList[PartPtr->_EOSID-1];
				
				if(EOSPtr->Type()==enIPCS)
				{
					((CSPHIPCS *)EOSPtr)->UpdateP(SPHPtPtr);
				}
			}
		}
	}
}
