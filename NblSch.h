/******************************************************************************
!
! Purpose:to get the _PtPairList
!
! Description:neighbour searching method
!
! Notes:
!
!******************************************************************************
!
!$Id: NBLSCH,v 2.0 (date)03.08 Richard LIU $
!
!Copyright: (c) 
!******************************************************************************/

#ifndef _NBLSCH_H_
#define _NBLSCH_H_

#include "Region.h"
#include "BasePt.h"
#include "PtPair.h"
#include "Box.h"
#include "PtMshPair.h"
#include "PtOperate.h"
#include <iostream>
#include <vector>
using namespace std;

class CNblSch
{
public:
	CNblSch(size_t CellNumx=1,size_t CellNumy=1);

	~CNblSch();

	void ResizeBox(size_t CellNumx,size_t CellNumy);

	void GetNbl(CRegion& Region);
	void GetMshNbl(CRegion& Region);

	void Clear(CRegion& Region);

	CBox _Box;//11.20
private:
  
	std::vector<CBasePt *> _PtNbList;

  void PushPtPair(vector<CPtPair> & PtPairList, CBasePt * PtiPtr, CBasePt * PtjPtr, enPTPAIRTYPE PtPairType, double xij, double yij, double r2);//used to push back particle pairs into PtPairList

}; 

#endif

/******************************************************************************
!
!
!******************************************************************************/
