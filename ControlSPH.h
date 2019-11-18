/******************************************************************************
!
! Purpose:to control all the variables of the processes
!
! Description:contain all of the control variables of the SPH particles
!
! Notes:
!
!******************************************************************************
!
!$Id: CONTROLSPH,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _CONTROLSPH_H_
#define _CONTROLSPH_H_

#include <string>
#include <iostream>
#include <sstream>
#include <iosfwd>

class CControlSPH
{
public:

	CControlSPH();

	~CControlSPH();

	unsigned int _RunMod;//RunMod=0 Implict; RumMod=1 Explicit

	unsigned int _StartStep;//0-Start from timeStep 0;esle-Start from StartStep

	unsigned int _VSL;

	unsigned int _SPHAV;

	unsigned int _SPHPV;//����ճ�Կ��Ʋ�����0��ţ��ճ�� 1��������ճ��

	unsigned int _SPHAS;

	unsigned int _SPHST;

	//unsigned int _SPHAD;

	unsigned int _SPHIPCS;


	unsigned int _PerdBnd;//0 ���������Ա߽磬1 ���������Ա߽磨����Poiseuille���ȣ�

	double _XSPHEpsilon;//XSPH�Ĳ���

	double _Cs;//����

	double _DeltaTCoeff;

	double _DeltaT;

	double _FinalTime;

	unsigned int _DensRenormSteps;//���ٲ��ܶ��ع�һ��

	unsigned int _OutputSteps;

	unsigned int _CellNumx;
	unsigned int _CellNumy;
	unsigned int _CellNumz;

	double _AVAlpha;
	double _AVBeta;
	double _AVEta;

	double _ASEpsilon1;
	double _ASEpsilon2;
	double _ASDeltaD;

	double _ADDelta;//Artificial diffusion (�˹���ɢ�A Colagrossi)

	double _Xcrmin;
	double _Xcrmax;
	double _Ycrmin;
	double _Ycrmax;
	double _Zcrmin;
	double _Zcrmax;


	//��������������ķ�Χ�������������Ѽ��㣩��ֻ��Y����
	double _STYmin;
	double _STYmax;


	//��������������������������Щ��

	unsigned int _CSPMIflag1;//CSPM�㷨�����˺�����ֵʱ�ķ�ĸ��1 �����Ѿ����㣬0 ��û����
	unsigned int _CSPMIflag2;//CSPM�㷨��������������ֵʱ��ϵ������

	//��ֵ�õ�������Ϣ
	double _MshMinX;
	double _MshMaxX;

	double _MshMinY;
	double _MshMaxY;

	double _MshSize;

	std::string _InfileName;//�����ļ�������


	double _PerdBndMinX;//�����Ա߽���������Сx������Poiseuille��
	double _PerdBndMaxX;

	unsigned int _IFlagRemesh;//�Ƿ������Remeh

	double _TDamp;//Dampʱ�䣬XYHu 2012Equ��13��

	double _PtShftBeta;//particle shift coefficient, beta

};

#endif

/******************************************************************************
!
!
!******************************************************************************/
