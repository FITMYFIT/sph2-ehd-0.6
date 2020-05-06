#include "SPHEqu.h"

CSPHEqu::CSPHEqu()
{
}

CSPHEqu::~CSPHEqu()
{
}


void CSPHEqu::Accelerate(CRegion &Region,double DeltaT,unsigned int TimeSteps)
{

  // if(Region._ControlSPH._RunMod==0)//隐式算法时，需要计算非稳态项
  //   {
  //     _UnsteadyTerm.Solve(Region);
  //   }

  // _SPHPressure.Solve(Region);

  // _SPHViscocity.Solve(Region,TimeSteps);

  // // _ExtForce.Solve2(Region,TimeSteps);//加入Damp过程的外力施加，XY Hu2012 Equ.(13)

  // _BndForce.Solve(Region);

  // if(Region._ControlSPH._SPHST==1)
  //   {
  //     _SPHSTForce.Solve(Region);
  //   }

  // if(Region._ControlSPH._SPHAS==1)
  //   {
  //     _SPHAStress.Solve(Region);
  //   }

  // if(Region._ControlSPH._SPHAV==1)
  //   {
  //     _SPHAVForce.Solve(Region);
  //   }
  // if(Region._ControlSPH._SPHAV==2)//Monaghan 1997型人工粘性
  //   {
  //     _SPHAVForce.Solve2(Region);
  //   }

  // if(Region._ControlSPH._VSL==1)
  //   {
  //     _SPHSmoothingEqu.Solve(Region);
  //   }

  //_SPHShearingForce.Solve(Region,DeltaT);

  //_SPHViscoelastic.Solve(Region,DeltaT);

  _SPHEHD.Solve4(Region, TimeSteps);//2019.09.23


  double RTC=1.0e-8;
  if(Region._ControlSPH._RunMod==0)
    {
      _EquSolve.BICGSolve(Region._SPHAm,Region._SPHu,Region._SPHbu,20,RTC);
      _EquSolve.BICGSolve(Region._SPHAm,Region._SPHv,Region._SPHbv,20,RTC);
    }
}
