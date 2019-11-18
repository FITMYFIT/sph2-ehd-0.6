/******************************************************************************
!
! Purpose: to define a type of particles
!
! Description:the type of SPH particles
!
! Notes:
!
!******************************************************************************
!
!$Id: SPHPT,v 2.0 (date)03.08 CHEN Fuzhen $
!
!Copyright: (c) 2010 by the Second Artillery Engineering College
!******************************************************************************/

#ifndef _SPHPT_H_
#define _SPHPT_H_

#include "BasePt.h"
#include <vector>

class CSPHPt : public CBasePt
{
public:

	CSPHPt(unsigned int ID,unsigned int PID);

	~CSPHPt();

	//the CSPM parameters
	double _Fxx;
	double _Fxy;
	double _Fxz;
	double _Fyx;
	double _Fyy;
	double _Fyz;
	double _Fzx;
	double _Fzy;
	double _Fzz;

	double _Cx;
	double _Cy;
	double _Cz;

	//the SCS parameters
	double _Gx;
	double _Gy;
	double _Gz;
	double _Gxx;
	double _Gxy;
	double _Gxz;
	double _Gyx;
	double _Gyy;
	double _Gyz;
	double _Gzx;
	double _Gzy;
	double _Gzz;
	double _Mm;
	double _Mmx;
	double _Mmy;
	double _Mmz;
	double _Mx;
	double _My;
	double _Mz;
	double _Mxx;
	double _Mxy;
	double _Mxz;
	double _Myy;
	double _Myx;
	double _Myz;
	double _Mzx;
	double _Mzy;
	double _Mzz;


	//the IPCS parameters
	double _deltap;
	double _deltap1;

	double _deltau;
	double _deltav;
	double _deltaw;

	double _sfx;
	double _sfy;
	double _sfz;

	//the shearing force parameters
	double _Rxy;
	double _Rxz;
	double _Ryz;
	double _J2;

	//the viscoelastic model parameters
	double _kxx;
	double _kxy;
	double _kxz;
	double _kyx;
	double _kyy;
	double _kyz;
	double _kzx;
	double _kzy;
	double _kzz;

	double _sigmaxx;
	double _sigmaxy;
	double _sigmaxz;
	double _sigmayy;
	double _sigmayz;
	double _sigmazz;

	//the parameters used by shearing force and viscoelastic model
	double _depsilonxx;
	double _depsilonxy;
	double _depsilonxz;
	double _depsilonyy;
	double _depsilonyz;
	double _depsilonzz;

	double _dsxx;
	double _dsxy;
	double _dsxz;
	double _dsyy;
	double _dsyz;
	double _dszz;
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
