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

#ifndef _MESH_H_
#define _MESH_H_

#include "Node.h"
#include "BasePt.h"
#include <vector>

class CMesh
{
public:

	CMesh();

	~CMesh();

	//enLOCALTYPE _LocalType;

	unsigned int _Localx;//x�����Ƿ��ڱ�Ե1����һ�㣻2����2�㣻3����3�㣻0�ڲ�
	unsigned int _Localy;



	unsigned int _Iflag;//Iflag=0 ���ڼ�����Iflag=1 �ڼ�������

	unsigned int _ID;

	unsigned int _PID;

	unsigned int _ID2;//�������������ӱ�ţ�ÿ����ͬ

	unsigned int _ID3;//������ż���������ӱ�ţ�ÿ����ͬ


	double _x;

	double _y;

	double _z;



	double _h;//smoothing length
	double _H;

	double _r;
	double _rr;

	double _Volume;

	double _C0;//��ʼɫֵ�����������ظ��ˣ����޸��˼�����ų����ȥ��
	double _C;//��ֵ���ɫֵ

	//the surface normal
	double _nx,_ny,_nz;
	double _Nnx,_Nny,_Nnz;	
	double _Curvature;

	CNode * _MeshNodList[4];//��ɸ������4���ڵ�

	unsigned int _NebrMeshID[4];//�������4����������������ID���������������У����ȱ�ĸ����򽫸�ID��0

	double _m;	//mass

	double _rho;//density

	double  _rho0;

	//pressure
	double _p;

	double _u;//x component
	double _v;//y component
	double _w;//z component

	//Ӧ��Remesh�㷨
	double _Mometu;//x����Ķ���
	double _Mometv;


	double _CSPMCoef;//CSPM����������ֵʱ�ķ�ĸ

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
