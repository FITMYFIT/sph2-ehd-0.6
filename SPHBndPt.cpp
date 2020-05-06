#include "SPHBndPt.h"

CSPHBndPt::CSPHBndPt(unsigned int ID,unsigned int PartID)
:CBasePt(enBndPt,ID,PartID),_nx(0.0),_ny(0.0),_nz(0.0),_area(0.0)
{
}

CSPHBndPt::~CSPHBndPt()
{
}
