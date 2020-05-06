/******************************************************************************
!
! Purpose:as the door of input and output
!
! Description:contain all the command of input and output
!
! Notes:
!
!******************************************************************************
!
!$Id: KFILE,v 2.0 (date)03.09 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _REMESH_H_
#define _REMESH_H_

#include <vector>
#include <cmath>
#include"Region.h"
#include"SPHPt.h"
#include"Part.h"
#include "NblSch.h"//11.20
#include "Matrix.h"
#include "Mesh.h"
#include "GetKnlList.h"

using namespace std;

class CRemesh
{
public:
	CRemesh();

	~CRemesh();

	void Solve(CRegion & Region);

private:

	CNblSch _GetMshNeigb;

	CGetKnlList _GetMshKnlList;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/

