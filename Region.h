/******************************************************************************
!
! Purpose:to offer a public region for all parameters
!
! Description:contain all parameters
!
! Notes:
!
!******************************************************************************
!
!$Id: REGION,v 2.0 RICHARD LIU $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _REGION_H_
#define _REGION_H_

#include <vector>
#include "BasePt.h"
#include "PtPair.h"
#include "Knl.h"
#include "Force.h"
#include "Part.h"
#include "Section.h"
#include "EOS.h"
//#include "Mat.h"
#include "ControlSPH.h"
#include "StatData.h"
#include "Matrix.h"
#include "Mesh.h"
#include "Node.h"
#include "PtMshPair.h"

class CRegion
{
public:
  CRegion();

  ~CRegion();

  std::vector<CBasePt> _PtList;

  std::vector<CMesh> _MeshList;

  std::vector<CNode> _NodeList;

  std::vector<CPtPair> _PtPairList;

  std::vector<CPtPair> _PtPairList2;

  std::vector<CPtMshPair> _PtMshPairList;//����-�����

  std::vector<CKnl> _KnlList;

  std::vector<CKnl> _KnlList2;

  std::vector<CKnl> _PtMshKnlList;//��Ӧ��Pt-Msh���Ӷ�

  std::vector<CForce *> _ExtForceList;//��λ�����

  std::vector<CPart> _PartList;

  std::vector<CSection *> _SectionList;

  std::vector<CEOS *> _EOSList;

  std::vector<unsigned int> _CalList;

  std::vector<CMatrix> _SPHAm;//��ʽ�ⷨ��ϵ������

  std::vector<double> _SPHbu;//��ʽ�ⷨ��Դ��
  std::vector<double> _SPHbv;
  std::vector<double> _SPHbw;

  std::vector<double> _SPHu;
  std::vector<double> _SPHv;
  std::vector<double> _SPHw;

  CControlSPH _ControlSPH;

  CStatData _StatDataList;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
