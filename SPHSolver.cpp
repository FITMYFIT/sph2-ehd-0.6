#include "SPHSolver.h"

CSPHSolver::CSPHSolver()
:_Ttime(0.0),_TimeSteps(0)
{
}

CSPHSolver::~CSPHSolver()
{
}

void CSPHSolver::Input()
{

	_Region._ControlSPH._InfileName=
		/*"static tank"*/
		// "Lid_driven";
    // "ehdplanercase2";
    // "EHDBulkRelax";
    //  "EHDIsoCondCylinder";
    "ehddrop";
  //"dambreak-ETSIN-3";
  //"dambreak-han";
  /*"shelldrop"*/
  //"movedrop";
  //"Lid-driven";
  //"Movesquare-s";
  //"mould";
  //"control-disk2";

	_KFile.Input(_Region);

  //_Model.EHDPlannar(_Region);//model of ehd plannar, lopez 2011 4.1
  // _Model.EHDBulkRelax(_Region);//model of ehd bulk relaxation, lopez 2011 4.2.1
  //_Model.EHDIsoCondCylinder(_Region);//model of ehd isolated conduction cylinder
   _Model.EHDDrop(_Region);//model of ehd droplet, lopezf 2011 4.3

   // _CalBndNorm.Solve(_Region);//������CSFģ�͵ķ�������߽����ӵķ���

	//_KFile.InputMesh(_Region);

	//_KFile.OutMsh(_Region,0,0);
  
	_KFile.OutTecplot2(_Region,0,_Region._ControlSPH._InfileName+"_Init");
  //_KFile.OutVTK(_Region,0,"ehdplanner-init");
  
  cout<<"-----------------------------------"<<endl;

  cout<<"The total particle number is"<<_Region._PtList.size()<<endl;

	_NblSch.ResizeBox(_Region._ControlSPH._CellNumx,_Region._ControlSPH._CellNumy);

	//_ExtGrid.Resize(2*_Region._ControlSPH._CellNumx,2*_Region._ControlSPH._CellNumy,2*_Region._ControlSPH._CellNumz);
}

void CSPHSolver::Output()
{
	_KFile.OutTecplot(_Region,_TimeSteps);
}

void CSPHSolver::Run()
{
  _TimeSteps=_Region._ControlSPH._StartStep;

  _Ttime=_Region._ControlSPH._DeltaT*_TimeSteps;

  for(;_Ttime<_Region._ControlSPH._FinalTime;)
	{
    _StartT=clock();
    
    cout<<"the time steps is"<<_TimeSteps<<endl;

    cout<<"the total time is"<<_Ttime<<endl;

    ////�����Ӳ�ֵ�������ϣ���ֻ��ֵһ����������Ч��
    //if (_TimeSteps==0)
    //{
    //	_Region._ControlSPH._IFlagRemesh=0;
    //	_CalculateRange._UpdateRange(_Region,_TimeSteps);
    //	_Remesh.Solve(_Region);
    //	_Region._ControlSPH._IFlagRemesh=1;
    //}

    cout<<"0.Preparation: Update Calculation Range & Initialization."<<endl;

    _CalculateRange._UpdateRange(_Region,_TimeSteps);

    _SPHInit.Solve(_Region,_TimeSteps);//��ʼ��
    
    cout<<"0.Preparation has been done."<<endl;

    cout<<"the number of particles involved in calculation is"<<_Region._CalList.size()<<endl;

    cout<<"1.start neighbor particles searching."<<endl;

    if(  _TimeSteps==_Region._ControlSPH._StartStep)//only for cases particle not move, particle pairs donot vary
      {

        _NblSch.GetNbl(_Region);     

        cout<<"2.start getting knllist."<<endl;

        _GetKnlList.GetKnlList(_Region);

        cout<<"2.knllist has been got."<<endl;
      }
    cout<<"1.neighbor searching has been done (total neighbor pairs:"<<_Region._PtPairList.size()<<")"<<endl;
    //��ʱ��
    //---------------------------------
    //_IdentPtLocal.IdentPtLocal(_Region);//��ʱ�ģ���һ������С���������Ľ��

    //-----------------------

    _DeltaT.GetDeltaT(_Region);

    cout<<"DeltaT:"<<_DeltaT._DeltaT1<<endl;

    cout<<"3.start calculating continuity equation."<<endl;

    //_ContinuityEqu.Solve(_Region,_TimeSteps);//�����Ҫ����������˹���ɢ�Ȳ���

    cout<<"3.continuity equation has been solved."<<endl;

    cout<<"4.start solving equation of station"<<endl;

    //_SPHEOS.Solve(_Region);

    cout<<"4.equation of station has been solved."<<endl;

    //��ֵ�õ�Dummy���ӵ�ѹ������������ܶ�
    //  _GetDumProperty.Solve(_Region);

    cout<<"5.start solving momentum equation."<<endl;

    _SPHEqu.Accelerate(_Region,_DeltaT._DeltaT1,_TimeSteps);

    cout<<"5.momentum equation has been solved."<<endl;

    cout<<"6.start update the particle position."<<endl;

    if(_Region._ControlSPH._RunMod==1)//��ʽ����ʱ��ҪLeapFrog
      {
        _UpdatePosition.LeapFrogUpdate(_Region,_DeltaT._DeltaT,_DeltaT._DeltaT1,_TimeSteps);
        }
      else//��ʽ�ƽ���������
        {
          _UpdatePosition.ImplicitUpdate(_Region);
        }
      cout<<"6.particle position has been updated."<<endl;

      _Ttime+=_DeltaT._DeltaT1;

      _TimeSteps++;

      cout<<_TimeSteps<<" time steps' calculation has been done"<<endl;

      _EndT=clock();

      cout<<"Run time is :"<<(double)(_EndT-_StartT)/CLOCKS_PER_SEC<<"s"<<endl;

      cout<<"----------------------------------------------------------"<<endl;

      if(_TimeSteps%(_Region._ControlSPH._OutputSteps)==0)
        {
          // Output();
        }

      //if particles not move, no need to clear particle pair list
      // _NblSch.Clear(_Region);

      //_GetKnlList.ClearKnlList(_Region);

      _CalculateRange._Clear(_Region);
	}
}

void CSPHSolver::Done()
{
	_NblSch.Clear(_Region);

	_GetKnlList.ClearKnlList(_Region);

	_KFile.Clear(_Region);

	if(_Region._ControlSPH._RunMod==0)//��ʽ���㣬��Ҫ��ռ��㷽�̵ľ���
	{
		_Region._SPHAm.clear();

		_Region._SPHu.clear();
		_Region._SPHv.clear();
		_Region._SPHw.clear();

		_Region._SPHbu.clear();
		_Region._SPHbv.clear();
		_Region._SPHbw.clear();
	}
}
