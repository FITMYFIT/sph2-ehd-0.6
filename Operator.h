/******************************************************************************
!
! Purpose:����BASIC��������������  
!
! Description:
!
! Notes:
!
!******************************************************************************
!
!$Id: operator,v 2.0 (date)03.18 $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _OPERTATOR_H_
#define _OPERTATOR_H_

#include "SPHPt.h"
#include "Region.h"
#include "Matrix.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <float.h>


using namespace std;

#define IsZero(a) (fabs(a) < 10.0 * DBL_MIN)
#define IsOne(a)  (fabs(a - 1.0) < 10.0 * DBL_EPSILON)

class COperator
{
public:

	COperator();

	~COperator();
	
	double SumSquare(const std::vector<double> &v);//1�����Ԫ�ص�ƽ���ͣ��ٿ���

	double SumFabs(const std::vector<double> &v);//2�����Ԫ�صľ���ֵ֮��
	
	double Multi_VV(const std::vector<double> &v1,const std::vector<double> &v2);//3.res=v1*v2
	
	void AddAsign_VV(std::vector<double> &v1,const std::vector<double> &v2);//4.v1+=v2
	
	void SubAsign_VV(std::vector<double> &v1,const std::vector<double> &v2);//5.v1-=v2

	void Subs_VV(const std::vector<double> &v1,const std::vector<double> &v2,std::vector<double> &vres);//6.res=v1-v2

	void Addi_VV(const std::vector<double> &v1,const std::vector<double> &v2,std::vector<double> &vres);//7.vres=v1+v2

	void Multi_SV(double S,const std::vector<double> &v,std::vector<double> &vres);//8.vres=S*v

	void Multi_MV(const std::vector<CMatrix> &A,const std::vector<double> &b,std::vector<double> &vres);//9.SRes=A*b

	void Multi_MtV(const std::vector<CMatrix> &A,const std::vector<double> &b,std::vector<double> &vres);//10.SRes=(A^T)*b

	void Reve2ndMat(double a11,double a12,double a21,double a22,
								double *ra11,double *ra12,double *ra21,double *ra22);//11.����׾���������
	void Reve3rdMat(double a11,double a12,double a13,double a21,double a22,double a23,double a31,double a32,double a33,
								double *ra11,double *ra12,double *ra13,double *ra21,double *ra22,double *ra23,double *ra31,double *ra32,double *ra33);//11.�����׾���������

	double MinEig2ndMat(double a11,double a12,double a21,double a22);//����׾������С����ֵ
	
	double LMax(double a1,double a2);//����������ֵ
	double LMin(double a1,double a2);//���������Сֵ

private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
