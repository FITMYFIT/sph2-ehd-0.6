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
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _NODE_H_
#define _NODE_H_

#include <vector>

enum enNODELOCALTYPE {enLeftBotm,enBotm,enRigtBotm,enRigt,enRigtUper,enUper,enLeftUper,enLeft,enIn};

class CNode
{
public:

	CNode();

	~CNode();

	enNODELOCALTYPE _LocalType;


	unsigned int _Iflag;//Iflag=0 ���ڼ�����Iflag=1 �ڼ�������

	unsigned int _ID;

	double _x;

	double _y;

	double _z;

	//the surface normal
	double _nx,_ny,_nz;
	double _Nnx,_Nny,_Nnz;	
	double _Curvature;

	double _VectMod;//����ģֵ

	std::vector<unsigned int> _MshIDList;//����ڵ����ڵ�����ı�ţ�ÿ���ڵ������4����������

	//double _m;	//mass

	//double _rho;//density

	//double  _rho0;

	////pressure
	//double _p;

	//double _u;//x component
	//double _v;//y component
	//double _w;//z component


	//double _h;//smoothing length
	//double _H;

	//double _r;
	//double _rr;

	//double _Volume;


	//double _uXSPHCoef;//XSPH�������������ƶ��ٶȵ�������
	//double _vXSPHCoef;
	//double _wXSPHCoef;

	//double _C0;//��ʼɫֵ�����������ظ��ˣ����޸��˼�����ų����ȥ��
	//double _C;//��ֵ���ɫֵ

	//double _CSPMCoef;//CSPM����������ֵʱ�ķ�ĸ

	//double _CSPMAxx;//CSPM���������ݶ�ʱ��ϵ������
	//double _CSPMAxy;
	//double _CSPMAxz;
	//double _CSPMAyx;
	//double _CSPMAyy;
	//double _CSPMAyz;
	//double _CSPMAzx;
	//double _CSPMAzy;
	//double _CSPMAzz;

	//double _CSPMminLambda;//CSPMϵ���������С����ֵ

private:
};



#endif

/******************************************************************************
!
!
!******************************************************************************/
