/******************************************************************************
!
! Purpose:consindering the shearingforce  
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHVISCOCITY,v 2.0 (date)03.18 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHVISCOSITY_H_
#define _SPHVISCOSITY_H_

#include "SPHPt.h"
#include "Region.h"
#include "SPHIPCS.h"
#include "WeaklyCompress.h"
#include "Matrix.h"
#include "GetKnlList.h"
#include "CSPM.h"
#include <cmath>

class CSPHViscosity
{
public:

	CSPHViscosity();

	~CSPHViscosity();

	void Solve(CRegion &Region,unsigned int TimeSteps);

private:

	std::vector<double> _Gradux;//三个方向的速度梯度，计算gamma的时候用
	std::vector<double> _Graduy;
	std::vector<double> _Graduz;
	std::vector<double> _Gradvx;
	std::vector<double> _Gradvy;
	std::vector<double> _Gradvz;
	std::vector<double> _Gradwx;
	std::vector<double> _Gradwy;
	std::vector<double> _Gradwz;

  CGetKnlList _GetKnlListV;

	CCSPM _CSPMV;

	void SolveEta(CRegion &Region,unsigned int TimeSteps);

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
