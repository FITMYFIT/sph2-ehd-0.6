/******************************************************************************
!
! Purpose:modelling part
!
! Description:model of the geometry
!
! Notes:
!
!******************************************************************************
!
!$Id: MODEL,v Richard LIU @DATA 20190926$
!
!Copyright: (c) GNU
!******************************************************************************/

#ifndef _MODEL_H_
#define _MODEL_H_

#include <vector>
#include "Region.h"
// #include "BasePt.h"
#include "SPHPt.h"
#include "common.h"

using namespace std;

class CModel
{
public:
  CModel();

  ~CModel();

  void EHDPlannar(CRegion & Region);//lopez 2011 4.1 problem

  void EHDBulkRelax(CRegion & Region);//Lopez 2011 4.2.1 problem

  void EHDIsoCondCylinder(CRegion & Region);//Lopez 2011 4.2.2 problem

  void EHDDrop(CRegion & Region);//lopez 2011 4.3 drop
  void EHDDrop2(CRegion & Region);//lopez 2011 4.3 drop , input model from truegrid model, set some variables here

private:

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
