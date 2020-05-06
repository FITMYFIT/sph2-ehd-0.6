/******************************************************************************
!
! Purpose:to get numerical values
!
! Description:SPH kernel solver
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHSOLVER,v  by Richard LIU $
!
!Copyright: (c)  by the Second Artillery Engineering University
!******************************************************************************/

#ifndef _SPHSOLVER_H_
#define _SPHSOLVER_H_

#include "Region.h"
#include "KFile.h"
#include "SPHEOS.h"
#include "NblSch.h"
#include "GetKnlList.h"
#include "SPHEqu.h"
#include "DeltaT.h"
#include "ControlSPH.h"
#include "SPHIPCSEqu.h"
#include "CalculateRange.h"
#include "ContinuityEqu.h"
#include "Operator.h"
#include "UpdatePosition.h"
#include "SPHInit.h"
#include "IdentPtLocal.h"
#include "Remesh.h"
#include "GetDumProperty.h"//2014.09.08
#include "CalBndNorm.h"//2014.09.09
#include "Model.h"//2019.09.26

#include <iostream>
#include <vector>
#include <ctime>
using namespace std;


class CSPHSolver
{
public:
  CSPHSolver();

  ~CSPHSolver();

  void Input();

  void Run();

  void Done();

private:

  double _Ttime;

  size_t _TimeSteps;

  void Output();

  CRegion _Region;

  CKFile _KFile;

  CDeltaT _DeltaT;

  CUpdatePosition _UpdatePosition;

  CSPHEOS _SPHEOS;

  CSPHIPCSEqu _SPHIPCSEqu;

  CContinuityEqu _ContinuityEqu;

  CSPHEqu _SPHEqu;

  CNblSch _NblSch;

  CGetKnlList _GetKnlList;

  //CBndCd _BndCd;

  CCalculateRange _CalculateRange;

  CSPHInit _SPHInit;//≥ı ºªØ

  CIdentPtLocal _IdentPtLocal;

  CRemesh _Remesh;

  CGetDumProperty _GetDumProperty;

  CCalBndNorm _CalBndNorm;

  CModel _Model;//directly modelling in the program

  clock_t _StartT,_EndT;

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
