#ifndef _INDENTPTLOCAL_H_
#define _INDENTPTLOCAL_H_

#include "Region.h"
#include "BasePt.h"
#include "Knl.h"
#include "PtPair.h"
#include "CSPM.h"
#include "Operator.h"
#include <vector>
#include <cmath>

class CIdentPtLocal
{
public:
	CIdentPtLocal();

	~CIdentPtLocal();

	void IdentPtLocal(CRegion & Region);

private:
	CCSPM _CSPMI;

	COperator _OperaterI;
};
#endif