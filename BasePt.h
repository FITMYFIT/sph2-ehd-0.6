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

  unsigned int _Iflag;//Iflag=0 ���ڼ�����Iflag=1 �ڼ�������

  unsigned int _ID;

  unsigned int _ID2;//�������������ӱ�ţ�ÿ����ͬ
  unsigned int _ID3;//������ż���������ӱ�ţ�ÿ����ͬ

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

  double _Cs;//���٣������ﶨ��һ����Ϊ���ڱ�����ʱʹ�ã���Region._ControlSPH����ʵ�Ѿ�������ˣ���KFile.c����Region._ControlSPH._Cs��ֵ��

  double _h;//smoothing length
  double _H;

  double _r;
  double _rr;

  double _Volume;

  //������������ǰ������SPHPt���еģ����ڶ���������
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

  double _VisEta;//���ճ��

  double _VisGamma;//��������


  double _uXSPHCoef;//XSPH�������������ƶ��ٶȵ�������
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


  double _C0;//��ʼɫֵ�����������ظ��ˣ����޸��˼�����ų����ȥ��
  double _C;//��ֵ���ɫֵ

  double _CSPMCoef;//CSPM����������ֵʱ�ķ�ĸ

  double _CSPMAxx;//CSPM���������ݶ�ʱ��ϵ������
  double _CSPMAxy;
  double _CSPMAxz;
  double _CSPMAyx;
  double _CSPMAyy;
  double _CSPMAyz;
  double _CSPMAzx;
  double _CSPMAzy;
  double _CSPMAzz;

  double _CSPMminLambda;//CSPMϵ���������С����ֵ


  //���²�����֪���Ǹ�ʲô�õģ��������˹�Ӧ����ʱ���õ��ˣ�������0
  double _sxx;
  double _sxy;
  double _sxz;
  double _syy;
  double _syz;
  double _szz;

  //�����ã����������ļ��ٶ�
  double _AccBndx;
  double _AccBndy;

  double _DistanceBnd;//���ڱ߽����ӣ��߽����Ӿ�߽�ķ�����룬����SMarrone�߽�ѹ��ȷ��

  double _uwall;//�����˶��ٶȣ����ڷǻ��Ʊ߽�
  double _vwall;//�����˶��ٶȣ����ڷǻ��Ʊ߽�

  //���²�������particle shift algorithm
  unsigned int _NumNegbor;//ÿ������֧�����ڵ�������,����particle shift
  double _r0;//��ֹ���룬����particle shift
  double _PtShftCoefu;
  double _PtShftCoefv;
  double _PtShftCoefw;

  //parameters used in EHD model 2019.09.12
  double _eEpsilon;//��, permitivity of a fluid, which is non-constant in interface regions of two fluid
  double _eKappa; //��, conductivity of a fluid, similar condition to ��
  double _ePhi; //��, electric potential
  double _eEx; //��, vector, electric field
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
