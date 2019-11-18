#include "Part.h"

CPart::CPart(unsigned int PID)
  :_PID(PID),_SECID(0),_MID(0),_EOSID(0),_PartPtList(0,NULL),_PartCalList(0,0)
{
}

CPart::~CPart()
{
}
