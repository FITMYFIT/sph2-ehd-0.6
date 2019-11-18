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
!$Id: BASEPT,v 2.0 (date)Richard LIU $
!
!Copyright: (c) GNU
!******************************************************************************/

#ifndef _BASEPT_H_
#define _BASEPT_H_

enum enPTTYPE {enSPHPt=1, enNULLPt,enBndPt,enGhost1Pt,enGhost2Pt,enDumPt,enEHDDumPt,enEHDBndPt};

class CBasePt
{
public:

  CBasePt();

  CBasePt(enPTTYPE pttype,unsigned int ID, unsigned int PID);

  CBasePt(unsigned int ID, unsigned int PID);

  ~CBasePt();

  enPTTYPE _Type;

  unsigned int _Iflag;//Iflag=0 不在计算域；Iflag=1 在计算域内

  unsigned int _ID;

  unsigned int _ID2;//进入计算域的粒子编号，每步不同
  unsigned int _ID3;//进入表张计算域的粒子编号，每步不同

  unsigned int _PID;

  double _x;

  double _y;

  double _z;


  double _m;	//mass

  double _rho;//density

  double  _rho0;

  double _mrho;//m/rho, since this expression is commonly used in SPH, we define this variable and calculate before the beginning of every time step, in SPHInit 2019.09.13

  //pressure
  double _p;

  double _u;//x component
  double _v;//y component
  double _w;//z component


  double _T;//temperature

  double _Cs;//声速（在这里定义一个是为了在变声速时使用，在Region._ControlSPH中其实已经定义过了，在KFile.c中用Region._ControlSPH._Cs赋值）

  double _h;//smoothing length
  double _H;

  double _r;
  double _rr;

  double _Volume;

  //各类增量，以前定义在SPHPt类中的，现在定义在这了
  double _drho;//d(rho)/d(t)
  double _du;
  double _dv;
  double _dw;

  double _rhoHf;//vaule of half steps
  double _uHf;	//velocity
  double _vHf;
  double _wHf;


  double _hHf;
  double _dh;

  double _VisEta;//表观粘度

  double _VisGamma;//剪切速率


  double _uXSPHCoef;//XSPH法求解出的粒子移动速度的修正量
  double _vXSPHCoef;
  double _wXSPHCoef;

  //energy
  double _e;
  double _de;
  double _eHf;

  //the surface normal
  double _nx,_ny,_nz;
  double _Nnx,_Nny,_Nnz;
  double _Curvature;
  double _midcurve;
  double _midcorrect;


  double _C0;//初始色值，定义在这重复了，待修改了计算表张程序后去掉
  double _C;//插值后的色值

  double _CSPMCoef;//CSPM修正函数插值时的分母

  double _CSPMAxx;//CSPM修正函数梯度时的系数矩阵
  double _CSPMAxy;
  double _CSPMAxz;
  double _CSPMAyx;
  double _CSPMAyy;
  double _CSPMAyz;
  double _CSPMAzx;
  double _CSPMAzy;
  double _CSPMAzz;

  double _CSPMminLambda;//CSPM系数矩阵的最小特征值


  //以下参数不知道是干什么用的，但计算人工应力的时候用到了，必须置0
  double _sxx;
  double _sxy;
  double _sxz;
  double _syy;
  double _syz;
  double _szz;

  //测试用，各力产生的加速度
  double _AccBndx;
  double _AccBndy;

  double _DistanceBnd;//用于边界粒子，边界粒子距边界的法向距离，用于SMarrone边界压力确定

  double _uwall;//壁面运动速度，用于非滑移边界
  double _vwall;//壁面运动速度，用于非滑移边界

  //以下参数用于particle shift algorithm
  unsigned int _NumNegbor;//每个粒子支持域内的粒子数,用于particle shift
  double _r0;//截止距离，用于particle shift
  double _PtShftCoefu;
  double _PtShftCoefv;
  double _PtShftCoefw;

  //parameters used in EHD model 2019.09.12
  double _eEpsilon;//ε, permitivity of a fluid, which is non-constant in interface regions of two fluid
  double _eKappa; //κ, conductivity of a fluid, similar condition to ε
  double _ePhi; //φ, electric potential
  double _eEx; //Ε, vector, electric field
  double _eEy;
  double _eRho;//rhoe,charge density
  double _deRho;//the increasement of the charge density, calculate from charge density continuity equ Lopez 2011 Equ(21)
  double _geEpsilonx;//gradient(epsilon),x direction
  double _geEpsilony;

  double _Fex,_Fey;//for debug, temp,2019.10.09

  private:
};

#endif

/******************************************************************************
!
!
!******************************************************************************/
