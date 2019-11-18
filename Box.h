/******************************************************************************
!
! Purpose:contribute to neighbour searching and to get the _PtNbList which is initial of searching
!
! Description:the box searching method
!
! Notes:
!
!******************************************************************************
!
!$Id: BOX,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _BOX_H_
#define _BOX_H_

#include <vector>
#include <iostream>

#include"Region.h"
#include"BasePt.h"

class CBox
{
public:
	CBox(size_t CellNumx=1,size_t CellNumy=1);

	~CBox();

	void Resize(size_t CellNumx,size_t CellNumy);
	
	void UpdateBox(CRegion & Region);
	
	void GetNbl(double x, double y, double r,std::vector<CBasePt *> & PtNbList);

	void GetNbl2(double x, double y, double r,std::vector<CBasePt *> & PtNbList,CRegion & Region);//周期性边界的粒子对搜索

	void GetNearBndPt(std::vector<CBasePt *> & BndPt);

private:
	void SetMPoint(CRegion & Region);
	
	void PutPt2Cell(CRegion & Region);
	
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

};

#endif

/******************************************************************************
!
!******************************************************************************/
