/******************************************************************************
!
! Purpose: 计算边界的法向，用CSF模型
!
! Description:the type of SPH particles
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHPT,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _CALBNDNORM_H_
#define _CALBNDNORM_H_

#include "BasePt.h"
#include "Region.h"
#include "NblSch.h"
#include "Operator.h"
#include "KFile.h"
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

class CCalBndNorm
{
public:

	CCalBndNorm();

	~CCalBndNorm();

	void Solve(CRegion & Region);
private:

	CRegion RegionBnd;//读入的辅助文件的数据都存储在这个Region中

	void InputBndAssit(CRegion & Region);

	void GetNblBnd(CRegion & Region);

	void CenterBnd(vector<double> &px,vector<double> &py,vector<double> &p);

	void PrimarySchBnd(double x, double y, double r,std::vector<CBasePt *> & PtNbList);

	void GetWBnd(double x,double y,double h,double distance2,double & W,double & Wx,double & Wy,double & Ww);

	void GetNormBnd(CRegion & Region);

	std::vector<CBasePt *> _PtNbList;

	COperator _OperatorBnd;


	std::vector<std::vector<CBasePt *> > _CellList;

	size_t _CellNumx;

	size_t _CellNumy;

	size_t _CellNumz;

	double _MaxPointx;

	double _MinPointx;

	double _MaxPointy;

	double _MinPointy;

	double _MaxPointz;

	double _MinPointz;

	double _MinPointr;

	double _CellSizex;

	double _CellSizey;

	double _CellSizez;

	CKFile _KFileBnd;

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
