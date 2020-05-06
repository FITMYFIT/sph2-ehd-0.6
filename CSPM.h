/******************************************************************************
!
! Purpose: calculate the cspm coefficient when used to calculate fun and derivative
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: CSPM ,v 1.0 (date)2019.09 Richard LIU $
!
!Copyright: (c) 2019 
!******************************************************************************/

#ifndef _CSPM_H_
#define _CSPM_H_

#include "BasePt.h"
#include "Region.h"
#include "Operator.h"
#include <vector>

class CCSPM
{
public:

	CCSPM();

	~CCSPM();

	void GetCSPMFunCorctCoef(CRegion &Region);//����CSPM����������ϵ��

	void GetCSPMGradCorctCoef(CRegion &Region);//����CSPM�Ժ�����������ʱ��ϵ������, inverse one

  void GetCSPMGradCorctCoef2(CRegion &Region);//ehd dummy and boundary  particle involed

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
