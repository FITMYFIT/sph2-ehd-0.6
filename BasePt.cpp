#include "BasePt.h"


CBasePt::CBasePt()
{
}

CBasePt::CBasePt(enPTTYPE pttype,unsigned int ID, unsigned int PID)
:_Type(pttype),_ID(ID),_PID(PID)
{
}


CBasePt::CBasePt(unsigned int ID, unsigned int PID)
:_ID(ID),_PID(PID)
{
}

CBasePt::~CBasePt()
{
}
