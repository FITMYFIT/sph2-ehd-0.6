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

	unsigned int _InvolvedFluidNum0;//��һʱ�䲽����������������
	unsigned int _InvolvedFluidNum;//����������������

	unsigned int _NullPtNum;

	unsigned int _BndPtNum;
	
	unsigned int _GhostPtNum;

	unsigned int _PtPairNum;

	unsigned int _PtPairNum2;//_PtPairList2���ӶԵ���Ŀ�����ڱ����������ʼ�������Ӷ�

	unsigned int _SPHPtPairNum;

	unsigned int _GhostPtPairNum;

	unsigned int _MshPtPairNum;

	unsigned int _KnlNum;

	unsigned int * _NumNeighbor;//ÿ�����ӵ��ٽ����ӵ���Ŀ

	unsigned int * _NeighborID;//ÿ�����ӵ��ٽ����ӵı��

	unsigned int _MaxNeighborNum;//ÿ��ʱ�䲽�����ӵ�����ٽ�������

	unsigned int _TotalNeighborNum;//ÿ�����ӵ��ٽ����Ӷ���Ŀ֮��

	//vector<unsigned int> _Index;//��ID2��ţ��ڲ������bpv�ȵı�ţ��������ڼ������ڲ�������������������֮��Ķ�Ӧ��ϵ

private:

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
