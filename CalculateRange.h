#ifndef _RANGEDRAW_H_
#define _RANGEDRAW_H_

#include"Region.h"
#include <iostream>

class CCalculateRange
{
public:

	CCalculateRange();

	~CCalculateRange();

	void _UpdateRange(CRegion& Region,unsigned int TimeSteps);

	void _Clear(CRegion& Region);

};

#endif