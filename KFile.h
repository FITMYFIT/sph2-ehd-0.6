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

#ifndef _KFILE_H_
#define _KFILE_H_

#include<vector>
#include"Region.h"
#include"SPHPt.h"
#include"Part.h"
#include "IdealGas.h"
#include "WeaklyCompress.h"
#include "SectionSPH.h"
#include "Force.h"
#include "SectionNULL.h"
#include "NblSch.h"//11.20
#include "SPHIPCS.h"
#include "Matrix.h"
#include "common.h"

using namespace std;

class CKFile
{
public:
  CKFile();

  ~CKFile();

  void Input(CRegion & Region);

  void InputMesh(CRegion & Region);//读插值用的背景网格函数

  void OutTecplot(CRegion & Region,unsigned int TimeSteps);
  
  void OutTecplot2(CRegion & Region, unsigned int TimeSteps, string ouputname);//same as OutTecplot, specify outputname

  void outTecplotEHDPLANNER(CRegion & Region, unsigned int TimeSteps, string outputname);//output for ehd planner test, Lopez 2011
  void outTecplotEHDBulkRelax(CRegion & Region, unsigned int TimeSteps, string outputname);//output for ehd bulk relaxation test, Lopez 2011
  void outTecplotIsoCondCylinder(CRegion & Region, unsigned int TimeSteps, string outputname);//output for ehd isolated conducting cylinder test, Lopez 2011
  void outTecplotEHDDrop(CRegion & Region, unsigned int TimeSteps, string outputname);//output for ehd droplet deformation test, Lopez 2011
  
  void OutTecplotMsh(CRegion & Region,unsigned int TimeSteps);

  void OutVTK(CRegion & Region, unsigned int TimeSteps, string ouputname);//output VTK file format, for paraview 2019.10.30
  
  void OutMsh(CRegion & Region,unsigned int TimeSteps,double CurTime );//按网格形式输出SPH背景网格

  void Clear(CRegion & Region);

  void OutVector(std::vector<double> &v);//输出矩阵与向量，测试矩阵方程计算是否正确用
  void OutMatrix(std::vector<CMatrix> &M,unsigned int dim);

  void OutPtNeighbour(CRegion & Region, CBasePt * PtPtr,unsigned int TimeSteps, string outputname);//for test, output particle pairs,2019.11.01

private:

  void Center(vector<double> & px,vector<double> & py,vector<double> & p); //六面体分解为五个四面体叠加求解体积与重心

  void InputRestart(CRegion & Region);

  double _InitSPH[10][4];//最多支持10个part,用来存模型的初始速度，密度和温度的，因为在重启程序中也要用，故定义在此
};

#endif

/******************************************************************************
!
!
!******************************************************************************/

