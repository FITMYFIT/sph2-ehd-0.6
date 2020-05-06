#include "PtPair.h"

CPtPair::CPtPair()
{
}


CPtPair::CPtPair(enPTPAIRTYPE Type, CBasePt * PtiPtr, CBasePt * PtjPtr,double driac2)
:_Type(Type),_PtiPtr(PtiPtr),_PtjPtr(PtjPtr),_driac2(driac2)
{
}

CPtPair::~CPtPair()
{
}
