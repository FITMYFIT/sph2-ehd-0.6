/******************************************************************************
!
! Purpose:to solve the N-S equations  
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHEQU,v 2.0 (date)03.09 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHEQU_H_
#define _SPHEQU_H_

#include "SPHPt.h"
#include "Region.h"
#include "SPHShearingForce.h"
#include "SPHViscoelastic.h"
#include "SPHViscosity.h"
#include "SPHPressure.h"
#include "ExtForce.h"
#include "SPHAStress.h"
#include "SPHAVForce.h"
#include "SPHSmoothingEqu.h"
#include "UnsteadyTerm.h"
#include "EquSolve.h"
#include "KFile.h"
#include "SPHSTForce.h"
#include "BndForce.h"
#include "SPHEHD.h" //2019.09.23

class CSPHEqu
{
public:
  CSPHEqu();

  ~CSPHEqu();

  void Accelerate(CRegion & Region,double DeltaT,unsigned int TimeSteps);

private:

  CSPHShearingForce _SPHShearingForce;

  //CSPHViscoelastic _SPHViscoelastic;

  CSPHViscosity _SPHViscocity;

  CSPHPressure _SPHPressure;

  CExtForce _ExtForce;

  CSPHSTForce _SPHSTForce;

  CSPHAStress _SPHAStress;

  CSPHAVForce _SPHAVForce;

  CSPHSmoothingEqu _SPHSmoothingEqu;

  CUnsteadyTerm _UnsteadyTerm;

  CEquSolve _EquSolve;

  CBndForce _BndForce;

  CSPHEHD _SPHEHD;

  CKFile _KFile;//临时的，测试BICG算法计算正确性的，输出矩阵，用后删除

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
