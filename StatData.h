/******************************************************************************
!
! Purpose:to extend various types of particles
!
! Description:contain the basic parameters of particles
!
! Notes:
!
!******************************************************************************
!
!$Id: BASEPT,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2014 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _STATDATA_H_
#define _STATDATA_H_

#include <vector>

using namespace std;

class CStatData
{
public:
	
	CStatData();
	
	~CStatData();
	
	unsigned int _PartNum;
	
	unsigned int _SPHNum;

	unsigned int _SPHFluidNum;

	unsigned int _InvolvedFluidNum0;//上一时间步进入计算域的粒子数
	unsigned int _InvolvedFluidNum;//进入计算域的粒子数

	unsigned int _NullPtNum;

	unsigned int _BndPtNum;
	
	unsigned int _GhostPtNum;

	unsigned int _PtPairNum;

	unsigned int _PtPairNum2;//_PtPairList2粒子对的数目，用于表面张力曲率计算的粒子对

	unsigned int _SPHPtPairNum;

	unsigned int _GhostPtPairNum;

	unsigned int _MshPtPairNum;

	unsigned int _KnlNum;

	unsigned int * _NumNeighbor;//每个粒子的临近粒子的数目

	unsigned int * _NeighborID;//每个粒子的临近粒子的编号

	unsigned int _MaxNeighborNum;//每个时间步的粒子的最大临近粒子数

	unsigned int _TotalNeighborNum;//每个粒子的临近粒子对数目之和

	//vector<unsigned int> _Index;//按ID2编号，内部存的是bpv等的编号，建立处于计算域内参与计算的粒子与总粒子之间的对应关系

private:

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
