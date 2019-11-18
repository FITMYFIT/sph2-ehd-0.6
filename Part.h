/******************************************************************************
!
! Purpose:to define the basic type of the part
!
! Description: each part is a collection of particles with same properties
!
! Notes:
!
!******************************************************************************
!
!$Id: PART,v  (date) Richard LIU $
!
!Copyright: (c)  
!******************************************************************************/

#ifndef _PART_H_
#define _PART_H_

#include <vector>
#include "SPHPt.h"

enum enPARTTYPE {enSPH=1, enNULL,enBnd,enGhost1,enGhost2,enDummy,enEHDDum,enEHDBnd};

class CPart
{
public:
  CPart(unsigned int PID=0);

  ~CPart();

  unsigned int _PID;

  unsigned int _SECID;

  enPARTTYPE _PartType;

  unsigned int _MID;

  unsigned int _EOSID;

  double _HdivDp;//�⻬���������Ӽ��ı�ֵ

  double _Mu;//ţ������ճ��ϵ����һ��part����һ������

  double _VisK;// ����ճ�ȣ���ţ������ʱ��VisK=Mu��VisN=1��������ʱ����Ӧ��k & n�����ϵ��������ָ����

  double _VisN;

  int _C0;//����ĳ�ʼ�궨ɫֵ

  double _STSigma;//��������ϵ��

  double _STHvsh;//���ڱ��ż���Ĺ⻬����������ֵ�ı�ֵ

  std::vector<CSPHPt *> _PartPtList;

  std::vector<unsigned int> _PartCalList;

  //parameters for ehd
  double _eEpsilon; //�� of a part
  double _eKappa;// �� of a part

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
