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

  double _HdivDp;//光滑长度与粒子间距的比值

  double _Mu;//牛顿流体粘度系数，一个part代表一种流体

  double _VisK;// 流体粘度，当牛顿流体时，VisK=Mu，VisN=1，幂律流时，对应于k & n（稠度系数和幂律指数）

  double _VisN;

  int _C0;//流体的初始标定色值

  double _STSigma;//表面张力系数

  double _STHvsh;//用于表张计算的光滑长度与正常值的比值

  std::vector<CSPHPt *> _PartPtList;

  std::vector<unsigned int> _PartCalList;

  //parameters for ehd
  double _eEpsilon; //ε of a part
  double _eKappa;// κ of a part

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
