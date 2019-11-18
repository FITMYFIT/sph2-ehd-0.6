/******************************************************************************
!
! Purpose:electrohydrodynamics model in  SPH
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHEHD,v 0.0 (date)2019.09.04 Richard LIU $
!
!Copyright: (c) 2019 by the Richard LIU
!******************************************************************************/

#ifndef _SPHEHD_H_
#define _SPHEHD_H_

#include "Region.h"
#include "BasePt.h"
#include "SPHPt.h"
#include "Part.h"
#include "common.h"
#include "CSPM.h"
#include "KFile.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

//using LASpack to calculate linear equation ,BICG or other
extern "C" {
#include <laspack/itersolv.h>
#include <laspack/rtc.h>
}
using namespace std;

class CSPHEHD
{
public:
	CSPHEHD();

	~CSPHEHD();

	void Solve(CRegion & Region, unsigned int Timesteps);// the scheme in Lopez 2011 equ 21 22, P10 procedure

  void Solve2(CRegion & Region, unsigned int Timesteps);// the scheme used in Basilisk,http://basilisk.dalembert.upmc.fr/src/ehd/implicit.h

  void Solve3(CRegion & Region, unsigned int Timesteps);// based on solve2, calculate drhe instead of rhoe n+1
private:

  void output(double ** a, double * b, int n);//for debug, ouput coeffecient matrix for nabla (nabla(epsilon phi))=rhoe

  void output2(QMatrix & a, Vector & x, Vector & b ,string outputname);//for debug, ouput coeffecient matrix for nabla (nabla(epsilon phi))=rhoe, laspack format

  void outputV(Vector & x, string outputname);//output vector, laspack type

  CCSPM _CSPMEHD;//for cspm corrective

  void InterpEHDDumPropterty(CRegion & Region,unsigned int TimeSteps);//interpolate phi of ehd dummy particle from ehd dummy and fluid particles

  void InterpEHDBndProperty(CRegion & Region,unsigned int TimeSteps);//interpolate E of ehd boundary particles from fluid particles
  // CSPHSTForce _SPHSTForce;

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
