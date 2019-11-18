#include "KFile.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
// #include <sys/stat.h>//used to make directry in linux
// #include <sys/types.h>


CKFile::CKFile()
{
}

CKFile::~CKFile()
{
}

void CKFile::Input(CRegion & Region)
{
	string keyword;
	string line;
	vector<double> X,Y,Z;
	float value;
	unsigned int PID,SECID,EOSID,ID;
	char tch;

	ifstream infile;
	ostringstream infilename;
	infilename<<Region._ControlSPH._InfileName<<".k";
	infile.open(infilename.str().data(),ios::in);
	if (!infile)
	{
		cout<<"KFile could not be open."<<endl;
		return;
	}

	else
	{
		cout<<"Start Reading KFile."<<endl;
	}


	while(!infile.eof())
	{
		infile>>keyword;

		if(keyword=="*RUN_MODE")
		{
			getline(infile,line);

			infile>>value;
			Region._ControlSPH._RunMod=size_t(value);//0-Implicit 1-Explicit
		}

		if(keyword=="*START_STEP")
		{
			getline(infile,line);

			infile>>value;
			Region._ControlSPH._StartStep=size_t(value);//0-从头开始算；othersteps-从某步开始计算
		}

		CPart PartI;
		unsigned int PtTypeID;
		if(keyword=="*PART")
		{
			getline(infile,line);

			for(;;)
			{
				infile>>tch;
				infile.seekg(-1,ios::cur);
				if((tch=='*')||(tch=='$'))	break;

				infile>>value;
				infile.ignore();
				PID=size_t(value);
				PartI._PID=PID;//1 PID

				infile>>value;
				infile.ignore();
				PtTypeID=size_t(value);//2

				infile>>value;
				infile.ignore();
				PartI._MID=size_t(value);//3

				infile>>value;
				infile.ignore();
				PartI._EOSID=size_t(value);//4

				infile>>value;
				infile.ignore();
				PartI._C0=int(value);//5初始色值，表张计算使用

				infile>>value;
				infile.ignore();
				PartI._HdivDp=double(value);//6

				infile>>value;
				infile.ignore();
				PartI._VisK=double(value);//7 k 值，牛顿流体时，k=Mu,n=1;幂律流时对应于稠度系数和幂指数

				infile>>value;
				infile.ignore();
				PartI._VisN=double(value);//8 n 值，牛顿流体时，k=Mu,n=1;幂律流时对应于稠度系数和幂指数

        infile>>value;
				infile.ignore();
        PartI._eKappa=double(value);//9 Κ, the conductivity of the part

        infile>>value;
        // infile.ignore();
        PartI._eEpsilon=double(value);//10 ε，the permitivity of the part


        //判断该Part属于哪种类型的粒子
        if(PtTypeID==1)
					PartI._PartType=enSPH;
				if(PtTypeID==2)
					PartI._PartType=enNULL;
				if(PtTypeID==3)
					PartI._PartType=enBnd;
				if(PtTypeID==4)
					PartI._PartType=enGhost1;
				if(PtTypeID==5)
					PartI._PartType=enGhost2;
				if(PtTypeID==6)
					PartI._PartType=enDummy;
        if(PtTypeID==7)
          PartI._PartType=enEHDBnd;
        if(PtTypeID==8)
          PartI._PartType=enEHDDum;

        Region._PartList.push_back(PartI);
			}
		}

		if(keyword=="*CONTROL_SPH1")
		{
			getline(infile,line);

			infile>>value;
			infile.ignore();
			Region._ControlSPH._SPHST=size_t(value);//1

			infile>>value;
			infile.ignore();
			Region._ControlSPH._SPHPV=size_t(value);//2

			infile>>value;
			infile.ignore();
			Region._ControlSPH._SPHAV=size_t(value);//3

			infile>>value;
			infile.ignore();
			Region._ControlSPH._SPHAS=size_t(value);//4

			infile>>value;
			infile.ignore();
			Region._ControlSPH._ADDelta=double(value);//5

			infile>>value;
			infile.ignore();
			Region._ControlSPH._XSPHEpsilon=double(value);//6

			infile>>value;
			infile.ignore();
			Region._ControlSPH._PerdBnd=size_t(value);//7

			infile>>value;
			infile.ignore();
			Region._ControlSPH._TDamp=double(value);//8. damp time ,X Y Hu 2012Equ 13

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._PtShftBeta=double(value);//9


			Region._ControlSPH._SPHIPCS=0;//人工将两个参数置0
			Region._ControlSPH._VSL=0;
		}

		if(keyword=="*CONTROL_SPH2")
		{
			getline(infile,line);

			infile>>value;
			infile.ignore();
			Region._ControlSPH._Cs=double(value);//1

			infile>>value;
			infile.ignore();
			Region._ControlSPH._DeltaT=double(value);//2

			infile>>value;
			infile.ignore();
			Region._ControlSPH._FinalTime=double(value);//3

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._OutputSteps=size_t(value);//4
		}

		if(keyword=="*CELL_NUM")
		{
			getline(infile,line);

			infile>>value;
			infile.ignore();
			Region._ControlSPH._CellNumx=size_t(value);//1

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._CellNumy=size_t(value);//2
		}

		if(keyword=="*CALCULATION_RANGE")
		{
			getline(infile,line);

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._Xcrmin=double(value);

			infile>>value;
			infile.ignore();
			Region._ControlSPH._Xcrmax=double(value);

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._Ycrmin=double(value);

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._Ycrmax=double(value);
		}

		if(keyword=="*SPH_AV")
		{
			getline(infile,line);

			if(Region._ControlSPH._SPHAV==1)
			{
				infile>>value;
				infile.ignore();
				Region._ControlSPH._AVAlpha=double(value);

				infile>>value;
				infile.ignore();
				Region._ControlSPH._AVBeta=double(value);

				infile>>value;
				//infile.ignore();
				Region._ControlSPH._AVEta=double(value);
			}
		}

		if(keyword=="*SPH_AS")
		{
			getline(infile,line);

			if(Region._ControlSPH._SPHAS==1)
			{
				infile>>value;
				infile.ignore();
				Region._ControlSPH._ASEpsilon1=double(value);

				infile>>value;
				infile.ignore();
				Region._ControlSPH._ASEpsilon2=double(value);

				infile>>value;
				//infile.ignore();
				Region._ControlSPH._ASDeltaD=double(value);
			}
		}

		//if(keyword=="*SPH_AD")
		//{
		//	getline(infile,line);

		//	if(Region._ControlSPH._SPHAD==1)
		//	{
		//		infile>>value;
		//		//infile.ignore();
		//		Region._ControlSPH._ADDelta=double(value);
		//	}
		//}

		if(keyword=="*DENSITY_RENORMSTEPS")
		{
			getline(infile,line);

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._DensRenormSteps=size_t(value);//0-No Renormlize；else-Renormlized Steps
		}

		if(keyword=="*PERIOD_BOUNDARY")
		{
			getline(infile,line);

			infile>>value;
			Region._ControlSPH._PerdBndMinX=double(value);//周期性边界条件的最小x，用于Poiseuille流计算

			infile>>value;
			Region._ControlSPH._PerdBndMaxX=double(value);//周期性边界条件的最小x，用于Poiseuille流计算
		}

		CIdealGas *IdealGasPtr;
		if(keyword=="*EOS_IDEAL_GAS")
		{
			getline(infile,line);

			double Cp,Mw,R;
			infile>>value;
			infile.ignore();
			EOSID=size_t(value);

			infile>>value;
			infile.ignore();
			Cp=double(value);

			infile>>value;
			infile.ignore();
			Mw=double(value);

			infile>>value;
			//infile.ignore();
			R=double(value);

			IdealGasPtr=new CIdealGas(Cp,Mw,R,EOSID);

			Region._EOSList.push_back(IdealGasPtr);
		}

		CWeaklyCompress *WeaklyCompressPtr;
		if(keyword=="*EOS_WEEKLY_COMPRESS")
		{
			getline(infile,line);

			double P0,gamma,ExtP;

			infile>>value;
			infile.ignore();
			EOSID=size_t(value);//1

			infile>>value;
			infile.ignore();
			P0=double(value);//2

			infile>>value;
			infile.ignore();
			gamma=double(value);//3

			infile>>value;
			ExtP=double(value);//4

			WeaklyCompressPtr=new CWeaklyCompress(P0,gamma,ExtP,EOSID);
			Region._EOSList.push_back(WeaklyCompressPtr);
		}

		if(keyword=="*SURFACE_TENSION")
		{
			getline(infile,line);

			if(Region._ControlSPH._SPHST==1)
			{
				for(unsigned int i=0;i!=Region._PartList.size();++i)
				{
					infile>>value;
					infile.ignore();
					PID=size_t(value);//1PID

					infile>>value;
					infile.ignore();
					Region._PartList[PID-1]._STSigma=double(value);//2sigma

					infile>>value;
					//infile.ignore();
					Region._PartList[PID-1]._STHvsh=double(value);//3STHvsh
				}
			}
		}

		CForce *ForcePtr;
		if(keyword=="*LOAD_BODYFORCE")
		{
			getline(infile,line);

			for(;;)
			{
				infile>>tch;
			  infile.seekg(-1,ios::cur);
			  if((tch=='*')||(tch=='$'))	break;

				infile>>value;
				infile.ignore();
				PID=size_t(value);
				ForcePtr=new CForce(PID);

				infile>>value;
				infile.ignore();
				ForcePtr->_fx=double(value);

				infile>>value;
				//infile.ignore();
				ForcePtr->_fy=double(value);


				Region._ExtForceList.push_back(ForcePtr);
			}
		}

		if(keyword=="*INITIAL_SPH")
		{
			getline(infile,line);

			for(unsigned int i=0;i!=Region._PartList.size();++i)
			{
				infile>>tch;
				infile.seekg(-1,ios::cur);
				if((tch=='*')||(tch=='$'))	break;

				infile>>value;
				infile.ignore();
				PID=size_t(value);//1PID

				infile>>value;
				infile.ignore();
				_InitSPH[PID-1][0]=double(value);//2u

				infile>>value;
				infile.ignore();
				_InitSPH[PID-1][1]=double(value);//3v

				infile>>value;
				infile.ignore();
				_InitSPH[PID-1][2]=double(value);//4rho

				infile>>value;
				infile.ignore();
				_InitSPH[PID-1][3]=double(value);//5T
			}
		}


		if(keyword=="*NODE")
		{
			if(Region._ControlSPH._StartStep==0)//从0时间步开始计算时读入该部分数据
			{
				for(;;)
				{
					infile>>tch;
					if((tch=='*')||(tch=='$'))
					{
						infile.seekg(-1,ios::cur);
						break;
					}
					while(tch!=',')
					{
						infile>>tch;
					}

					infile>>value;
					infile.ignore();
					X.push_back(value);

					infile>>value;
					infile.ignore();
					Y.push_back(value);

					for(int i=0;i<3;i++)
					{
						infile>>value;
						infile.ignore();
					}
				}
			}
		}

		if(keyword=="*ELEMENT_SHELL_THICKNESS")
		{
			if(Region._ControlSPH._StartStep==0)//从0时间步开始计算时读入该部分数据
			{

				for(;;)
				{
					vector<double>tPosx(4);
					vector<double>tPosy(4);
					vector<double>tPartile(3);

					CBasePt BasePt;

					infile>>tch;
					infile.seekg(-1,ios::cur);
					if((tch=='*')||(tch=='$'))	break;

					infile>>value;
					infile.ignore();
					ID=size_t(value);
					BasePt._ID=ID;

					infile>>value;
					infile.ignore();
					PID=size_t(value);
					BasePt._PID=PID;


					infile>>value;
					infile.ignore();
					ID=size_t(value);
					tPosx[0]=X[ID-1];
					tPosy[0]=Y[ID-1];

					infile>>value;
					infile.ignore();
					ID=size_t(value);
					tPosx[1]=X[ID-1];
					tPosy[1]=Y[ID-1];

					infile>>value;
					infile.ignore();
					ID=size_t(value);
					tPosx[2]=X[ID-1];
					tPosy[2]=Y[ID-1];

					infile>>value;
					infile.ignore();
					ID=size_t(value);
					tPosx[3]=X[ID-1];
					tPosy[3]=Y[ID-1];

					for(int i=0;i<4;i++)
					{
						infile>>value;
						infile.ignore();
					}

					Center(tPosx,tPosy,tPartile);

					BasePt._x=tPartile[0];
					BasePt._y=tPartile[1];
					BasePt._Volume=tPartile[2];

					BasePt._u=_InitSPH[PID-1][0];
					BasePt._v=_InitSPH[PID-1][1];

					BasePt._uwall=_InitSPH[PID-1][0];//暂时的，用于给壁面part的速度赋值，在GetBndProperty中用到了
					BasePt._vwall=_InitSPH[PID-1][1];//暂时的，用于给壁面part的速度赋值，在GetBndProperty中用到了

					BasePt._rho0=_InitSPH[PID-1][2];
					BasePt._rho =_InitSPH[PID-1][2];
					BasePt._T =_InitSPH[PID-1][3];

					BasePt._m=BasePt._rho0*BasePt._Volume;

					BasePt._h=Region._PartList[PID-1]._HdivDp*pow(BasePt._Volume,0.5);
					BasePt._r=2.0*BasePt._h;
					BasePt._C0=Region._PartList[PID-1]._C0;

					BasePt._Cs=Region._ControlSPH._Cs;//为每个粒子的声速赋初值

					if(Region._ControlSPH._SPHST==1)
					{
						BasePt._H=Region._PartList[PID-1]._STHvsh*BasePt._h;
						BasePt._rr=2*BasePt._H;
					}

					BasePt._p=0;
					BasePt._Iflag=1;
					switch(Region._PartList[PID-1]._PartType)
					{
				 	  case enSPH:
								BasePt._Type=enSPHPt;
								break;
				 	  case enNULL:
								BasePt._Type=enNULLPt;
								break;
				 	  case enBnd:
								BasePt._Type=enBndPt;
								break;
				 	  case enGhost1:
								BasePt._Type=enGhost1Pt;
								break;
						case enGhost2:
							BasePt._Type=enGhost2Pt;
							break;
						case enDummy:
							BasePt._Type=enDumPt;
							break;
						default:
								break;
					}

					//Region._PartList[PID-1]._PartPtList.push_back(&BasePt);

					Region._PtList.push_back(BasePt);
					tPosx.clear();
					tPosy.clear();
					tPartile.clear();
				}
			}
		}

		if(keyword=="*END")
		{
			break;
		}
	}


	if(Region._ControlSPH._StartStep!=0)//重启程序，从输出文件中读入模型数据
	{
		InputRestart(Region);
	}
}

void CKFile::Center(vector<double> &px,vector<double> &py,vector<double> &p)
{
	double x5,x6,y5,y6,s12,s13,s14,s23,s34,s1,s2,A1,A2,A,xc1,xc2,yc1,yc2;
	x5=(px[2]+px[3])/2;
	y5=(py[2]+py[3])/2;
	x6=(px[0]+px[2])/2;
  y6=(py[0]+py[2])/2;
	s12=sqrt((px[0]-px[1])*(px[0]-px[1])+(py[0]-py[1])*(py[0]-py[1]));
	s13=sqrt((px[0]-px[2])*(px[0]-px[2])+(py[0]-py[2])*(py[0]-py[2]));
	s14=sqrt((px[0]-px[3])*(px[0]-px[3])+(py[0]-py[3])*(py[0]-py[3]));
	s23=sqrt((px[1]-px[2])*(px[1]-px[2])+(py[1]-py[2])*(py[1]-py[2]));
	s34=sqrt((px[2]-px[3])*(px[2]-px[3])+(py[2]-py[3])*(py[2]-py[3]));
	s1=(s14+s13+s34)/2;
	s2=(s13+s12+s23)/2;
	A1=sqrt(s1*(s1-s14)*(s1-s13)*(s1-s34));
	A2=sqrt(s2*(s2-s12)*(s2-s23)*(s2-s13));
	xc1=(2*x5+px[0])/3;
	yc1=(2*y5+py[0])/3;
	xc2=(2*x6+px[1])/3;
	yc2=(2*y6+py[1])/3;
	A=A1+A2;
	p[0]=(xc1*A1+xc2*A2)/A;//输出重心坐标x
	p[1]=(yc1*A1+yc2*A2)/A;//输出重心坐标y
	p[2]=A;//输出面积
}

void CKFile::OutTecplot(CRegion & Region,unsigned int TimeSteps)
{
  ofstream outfile;
  ostringstream outfilename;
  outfilename<<"ZTEC-"<<Region._ControlSPH._InfileName<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".dat";
  outfile.open(outfilename.str().data(),ios::out);

  if (!outfile)
      {
        cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
        return;
      }


    outfile
    <<"VARIABLES= "
    <<"\"AXISX\" "
    <<"\"AXISY\" "
    //<<"\"AXISZ\" "
    <<"\"PID\" "
    <<"\"ID\" "
    <<"\"ID2\" "
    <<"\"DENSITY\" "
    //<<"\"DENSITY0\" "
    <<"\"VELOCITYX\" "
    <<"\"VELOCITYY\" "
    //<<"\"VELOCITYZ\" "
    <<"\"PRESSURE\" "
    <<"\"SMOOTHLENGTH\" "
    <<"\"PTSIZE\" "//particle size, for visualization simpilicity in tecplot
    <<"\"IFLAG\" "
    <<"\"CI\" "
    <<"\"EEPSILON\" "
    <<"\"EKAPPA\" "
    <<"\"ERHOE\" "
    <<"\"EPHI\" "
    <<"\"EEX\" "
    <<"\"EEY\" "
    <<endl;
  outfile<<"ZONE T="<<"\""<<Region._ControlSPH._InfileName<<" step "<<TimeSteps<<"\" I= "<<Region._PtList.size()<<", F=POINT"<<endl;//I= 后面必须有个空格，否则重启读数不行

  CSPHPt *SPHPtPtr;
  for (size_t iii=0;iii<Region._PtList.size();iii++)
    {
      SPHPtPtr=(CSPHPt *)(&Region._PtList[iii]);
      {
        outfile
          <<setiosflags(ios_base::scientific)
          <<setprecision(16)
          <<setw(26)<<SPHPtPtr->_x<<" "
          <<setw(26)<<SPHPtPtr->_y<<" "
          //<<setw(26)<<SPHPtPtr->_z<<" "
          <<setw(26)<<SPHPtPtr->_PID<<" "
          <<setw(26)<<SPHPtPtr->_ID<<" "
          <<setw(26)<<SPHPtPtr->_ID2<<" "
          <<setw(26)<<SPHPtPtr->_rho<<" "
          //<<setw(26)<<SPHPtPtr->_rho0<<" "
          <<setw(26)<<SPHPtPtr->_u<<" "
          <<setw(26)<<SPHPtPtr->_v<<" "
          //<<setw(26)<<SPHPtPtr->_w<<" "
          <<setw(26)<<SPHPtPtr->_p<<" "
          <<setw(26)<<SPHPtPtr->_h<<" "
          <<setw(26)<<SPHPtPtr->_h/(Region._PartList[SPHPtPtr->_PID-1]._HdivDp)<<" "
          <<setw(26)<<SPHPtPtr->_Iflag<<" "
          <<setw(26)<<SPHPtPtr->_C<<" "
          <<setw(26)<<SPHPtPtr->_eEpsilon<<" "
          <<setw(26)<<SPHPtPtr->_eKappa<<" "
          <<setw(26)<<SPHPtPtr->_eRho<<" "
          <<setw(26)<<SPHPtPtr->_ePhi<<" "
          <<setw(26)<<SPHPtPtr->_eEx<<" "
          <<setw(26)<<SPHPtPtr->_eEy<<" "

          <<endl;
      }
    }


  cout<<"Output file "<<"\"" <<outfilename.str() <<"\""<<" has been written."<<endl;
   
  outfile.close();
}
void CKFile::OutTecplot2(CRegion & Region,unsigned int TimeSteps, string outputname)
{
  ofstream outfile;
  ostringstream outfilename;
  outfilename<<outputname<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".dat";
  outfile.open(outfilename.str().data(),ios::out);

  if (!outfile)
    {
      cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
      return;
    }


  outfile
    <<"VARIABLES= "
    <<"\"AXISX\" "
    <<"\"AXISY\" "
    //<<"\"AXISZ\" "
    <<"\"PID\" "
    <<"\"ID\" "
    <<"\"ID2\" "
    <<"\"DENSITY\" "
    //<<"\"DENSITY0\" "
    <<"\"VELOCITYX\" "
    <<"\"VELOCITYY\" "
    //<<"\"VELOCITYZ\" "
    <<"\"PRESSURE\" "
    <<"\"SMOOTHLENGTH\" "
    <<"\"PTSIZE\" "//particle size, for visualization simpilicity in tecplot
    <<"\"IFLAG\" "
    <<"\"CI\" "
    <<"\"EEPSILON\" "
    <<"\"EKAPPA\" "
    <<"\"ERHOE\" "
    <<"\"EPHI\" "
    <<"\"EEX\" "
    <<"\"EEY\" "
    <<endl;
  outfile<<"ZONE T="<<"\""<<Region._ControlSPH._InfileName<<" step "<<TimeSteps<<"\" I= "<<Region._PtList.size()<<", F=POINT"<<endl;//I= 后面必须有个空格，否则重启读数不行

  double beta=3;//used for ehdplanner, lopez 2011 4.1
  CSPHPt *SPHPtPtr;
  for (size_t iii=0;iii<Region._PtList.size();iii++)
    {
      SPHPtPtr=(CSPHPt *)(&Region._PtList[iii]);
      {
        outfile
          <<setiosflags(ios_base::scientific)
          <<setprecision(16)
          <<setw(26)<<SPHPtPtr->_x<<" "
          <<setw(26)<<SPHPtPtr->_y<<" "
          //<<setw(26)<<SPHPtPtr->_z<<" "
          <<setw(26)<<SPHPtPtr->_PID<<" "
          <<setw(26)<<SPHPtPtr->_ID<<" "
          <<setw(26)<<SPHPtPtr->_ID2<<" "
          <<setw(26)<<SPHPtPtr->_rho<<" "
          //<<setw(26)<<SPHPtPtr->_rho0<<" "
          <<setw(26)<<SPHPtPtr->_u<<" "
          <<setw(26)<<SPHPtPtr->_v<<" "
          //<<setw(26)<<SPHPtPtr->_w<<" "
          <<setw(26)<<SPHPtPtr->_p<<" "
          <<setw(26)<<SPHPtPtr->_h<<" "
          <<setw(26)<<SPHPtPtr->_h/(Region._PartList[SPHPtPtr->_PID-1]._HdivDp)<<" "
          <<setw(26)<<SPHPtPtr->_Iflag<<" "
          <<setw(26)<<SPHPtPtr->_C<<" "
          <<setw(26)<<SPHPtPtr->_eEpsilon<<" "
          <<setw(26)<<SPHPtPtr->_eKappa<<" "
          <<setw(26)<<SPHPtPtr->_eRho<<" "
          <<setw(26)<<SPHPtPtr->_ePhi<<" "
          <<setw(26)<<SPHPtPtr->_eEx<<" "
          <<setw(26)<<SPHPtPtr->_eEy<<" "
          <<endl;
      }
    }


  cout<<"Output file "<<"\"" <<outfilename.str() <<"\""<<" has been written."<<endl;

  outfile.close();
}

void CKFile::outTecplotEHDBulkRelax(CRegion & Region,unsigned int TimeSteps, string outputname)
{
  ofstream outfile;
  ostringstream outfilename;
  outfilename<<outputname<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".dat";
  outfile.open(outfilename.str().data(),ios::out);

  if (!outfile)
    {
      cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
      return;
    }


  outfile
    <<"VARIABLES= "
    <<"\"AXISX\" "
    <<"\"AXISY\" "
    //<<"\"AXISZ\" "
    <<"\"PID\" "
    <<"\"ID\" "
    <<"\"ID2\" "
    <<"\"DENSITY\" "
    //<<"\"DENSITY0\" "
    <<"\"VELOCITYX\" "
    <<"\"VELOCITYY\" "
    //<<"\"VELOCITYZ\" "
    <<"\"PRESSURE\" "
    <<"\"SMOOTHLENGTH\" "
    <<"\"PTSIZE\" "//particle size, for visualization simpilicity in tecplot
    <<"\"IFLAG\" "
    <<"\"DRHOE\" "
    <<"\"Fex\" "
    <<"\"Fey\" "
    <<"\"EPSILON\" "
    <<"\"KAPPA\" "
    <<"\"RHOE\" "
    <<"\"RHOEExt\" "
    <<"\"PHI\" "
    <<"\"EX\" "
    <<"\"EY\" "
    // <<"\"PHIext\" "
    // <<"\"Eext\" "
    // <<"\"ErrorPHI\" "
    // <<"\"ErrorE\" "
    <<"\"NUMNEIGHBOR\" "
      
    <<endl;
  outfile<<"ZONE T="<<"\""<<Region._ControlSPH._InfileName<<" step "<<TimeSteps<<"\" I= "<<Region._PtList.size()<<", F=POINT"<<endl;//I= 后面必须有个空格，否则重启读数不行
  CSPHPt *SPHPtPtr;

   //calculate the exact values of electric field and electric potential, and the errors between simulation and exact
  double beta,eta;//

  vector <double> RhoeExt;
  vector<double> ErrorRhoe;
  RhoeExt.resize(Region._PtList.size(),0.0);
  ErrorRhoe.resize(Region._PtList.size(),0.0);
 
  double t=(TimeSteps)*Region._ControlSPH._DeltaT;
  double r2;
  double rhoe0;
  double a=0.05;
  double epsilon;
  double kappa;
  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      SPHPtPtr=(CSPHPt *)(&Region._PtList[i]);

      if(SPHPtPtr->_Type==enSPHPt)
        {
          r2=pow(SPHPtPtr->_x,2)+pow(SPHPtPtr->_y,2);
          rhoe0=pow(e,-r2/(2*a*a))/(a*sqrt(2*PI));

          kappa=Region._PartList[SPHPtPtr->_PID-1]._eKappa;
          epsilon=Region._PartList[SPHPtPtr->_PID-1]._eEpsilon;

          RhoeExt[i]=rhoe0*pow(e,-kappa*t/epsilon);
        }
    }
  

  for (size_t iii=0;iii<Region._PtList.size();iii++)
    {
      SPHPtPtr=(CSPHPt *)(&Region._PtList[iii]);
      {
        outfile
          <<setiosflags(ios_base::scientific)
          <<setprecision(16)
          <<setw(26)<<SPHPtPtr->_x<<" "
          <<setw(26)<<SPHPtPtr->_y<<" "
          //<<setw(26)<<SPHPtPtr->_z<<" "
          <<setw(26)<<SPHPtPtr->_PID<<" "
          <<setw(26)<<SPHPtPtr->_ID<<" "
          <<setw(26)<<SPHPtPtr->_ID2<<" "
          <<setw(26)<<SPHPtPtr->_rho<<" "
          //<<setw(26)<<SPHPtPtr->_rho0<<" "
          <<setw(26)<<SPHPtPtr->_u<<" "
          <<setw(26)<<SPHPtPtr->_v<<" "
          //<<setw(26)<<SPHPtPtr->_w<<" "
          <<setw(26)<<SPHPtPtr->_p<<" "
          <<setw(26)<<SPHPtPtr->_h<<" "
          <<setw(26)<<SPHPtPtr->_h/(Region._PartList[SPHPtPtr->_PID-1]._HdivDp)<<" "
          <<setw(26)<<SPHPtPtr->_Iflag<<" "
          <<setw(26)<<SPHPtPtr->_deRho<<" "
          <<setw(26)<<SPHPtPtr->_Fex<<" "
          <<setw(26)<<SPHPtPtr->_Fey<<" "
          <<setw(26)<<SPHPtPtr->_eEpsilon<<" "
          <<setw(26)<<SPHPtPtr->_eKappa<<" "
          <<setw(26)<<SPHPtPtr->_eRho<<" "
          <<setw(26)<<RhoeExt[iii]<<" "//exact phi, used for ehdplanner test, Lopez 2011 Table1
          <<setw(26)<<SPHPtPtr->_ePhi<<" "
          <<setw(26)<<SPHPtPtr->_eEx<<" "
          <<setw(26)<<SPHPtPtr->_eEy<<" "
           <<setw(26)<<SPHPtPtr->_NumNegbor<<" "//error e, used for ehdplanner test
          <<endl;
      }
    }

  //output max rhoe
  outfile.close();
  unsigned int IDCenter=2245;//ID OF THE CENTER PARTICLE, its has the max rhoe
  outfilename.str("");
  outfilename<<"MaxRhoe-BulkRelax"<<outputname<<".dat";
  outfile.open(outfilename.str().data(),ios::out|ios::app);
  if (!outfile)
    {
      cerr<<"Out put file centerid could not be open."<<endl<<"tecplot.dat"<<endl;
      return;
    }
  outfile
    <<setiosflags(ios_base::scientific)
    <<setprecision(16)
    <<setw(26)<<t<<" "
    <<setw(26)<<Region._PtList[IDCenter-1]._eRho<<" "
    <<setw(26)<<RhoeExt[IDCenter-1]<<" "
    <<endl;

  //output rhe time evolution of line y=0
  outfile.close();
  outfilename.str("");
  outfilename<<"y=0-Rhoe-BulkRelax"<<outputname<<TimeSteps<<".dat";
  outfile.open(outfilename.str().data(),ios::out);
  if (!outfile)
    {
      cerr<<"Out put file y=0 line could not be open."<<endl<<"tecplot.dat"<<endl;
      return;
    }

  
  outfile
    <<"VARIABLES= "
    <<"\"AXISX\" "
    <<"\"Rhoe\" "
    <<"\"RhoeExt\" "
    <<endl;
  outfile<<"ZONE T="<<"\""<<Region._ControlSPH._InfileName<<" t= "<<t<<"\" I= "<<67<<", F=POINT"<<endl;
  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      if(ISZERO(Region._PtList[i]._y))
        {
          outfile
            <<setiosflags(ios_base::scientific)
            <<setprecision(16)
            <<setw(26)<<Region._PtList[i]._x<<" "
            <<setw(26)<<Region._PtList[i]._eRho<<" "
            <<setw(26)<<RhoeExt[i]<<" "
            <<endl;
        }
    }

  

  RhoeExt.clear();

  cout<<"Output file "<<"\"" <<outfilename.str() <<"\""<<" has been written."<<endl;

  outfile.close();
}
void CKFile::outTecplotIsoCondCylinder(CRegion & Region,unsigned int TimeSteps, string outputname)
{
  ofstream outfile;
  ostringstream outfilename;
  outfilename<<outputname<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".dat";
  outfile.open(outfilename.str().data(),ios::out);

  if (!outfile)
    {
      cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
      return;
    }


  outfile
    <<"VARIABLES= "
    <<"\"AXISX\" "
    <<"\"AXISY\" "
    //<<"\"AXISZ\" "
    <<"\"PID\" "
    <<"\"ID\" "
    <<"\"ID2\" "
    <<"\"DENSITY\" "
    //<<"\"DENSITY0\" "
    <<"\"VELOCITYX\" "
    <<"\"VELOCITYY\" "
    //<<"\"VELOCITYZ\" "
    <<"\"PRESSURE\" "
    <<"\"SMOOTHLENGTH\" "
    <<"\"PTSIZE\" "//particle size, for visualization simpilicity in tecplot
    <<"\"IFLAG\" "
    <<"\"DRHOE\" "
    <<"\"Fex\" "
    <<"\"Fey\" "
    <<"\"EPSILON\" "
    <<"\"KAPPA\" "
    <<"\"RHOE\" "
    //  <<"\"RHOEExt\" "
    <<"\"PHI\" "
    <<"\"EX\" "
    <<"\"EY\" "
    // <<"\"PHIext\" "
    // <<"\"Eext\" "
    // <<"\"ErrorPHI\" "
    // <<"\"ErrorE\" "
    <<"\"NUMNEIGHBOR\" "
      
    <<endl;
  outfile<<"ZONE T="<<"\""<<Region._ControlSPH._InfileName<<" step "<<TimeSteps<<"\" I= "<<Region._PtList.size()<<", F=POINT"<<endl;//I= 后面必须有个空格，否则重启读数不行
  CSPHPt *SPHPtPtr;

  //calculate the exact values of electric field and electric potential, and the errors between simulation and exact
  double beta,eta;//

  vector <double> RhoeExt;
  vector<double> ErrorRhoe;
  RhoeExt.resize(Region._PtList.size(),0.0);
  ErrorRhoe.resize(Region._PtList.size(),0.0);
 
  double t=(TimeSteps)*Region._ControlSPH._DeltaT;
  for (size_t iii=0;iii<Region._PtList.size();iii++)
    {
      SPHPtPtr=(CSPHPt *)(&Region._PtList[iii]);
      {
        outfile
          <<setiosflags(ios_base::scientific)
          <<setprecision(16)
          <<setw(26)<<SPHPtPtr->_x<<" "
          <<setw(26)<<SPHPtPtr->_y<<" "
          //<<setw(26)<<SPHPtPtr->_z<<" "
          <<setw(26)<<SPHPtPtr->_PID<<" "
          <<setw(26)<<SPHPtPtr->_ID<<" "
          <<setw(26)<<SPHPtPtr->_ID2<<" "
          <<setw(26)<<SPHPtPtr->_rho<<" "
          //<<setw(26)<<SPHPtPtr->_rho0<<" "
          <<setw(26)<<SPHPtPtr->_u<<" "
          <<setw(26)<<SPHPtPtr->_v<<" "
          //<<setw(26)<<SPHPtPtr->_w<<" "
          <<setw(26)<<SPHPtPtr->_p<<" "
          <<setw(26)<<SPHPtPtr->_h<<" "
          <<setw(26)<<SPHPtPtr->_h/(Region._PartList[SPHPtPtr->_PID-1]._HdivDp)<<" "
          <<setw(26)<<SPHPtPtr->_Iflag<<" "
          <<setw(26)<<SPHPtPtr->_deRho<<" "
          <<setw(26)<<SPHPtPtr->_Fex<<" "
          <<setw(26)<<SPHPtPtr->_Fey<<" "
          <<setw(26)<<SPHPtPtr->_eEpsilon<<" "
          <<setw(26)<<SPHPtPtr->_eKappa<<" "
          <<setw(26)<<SPHPtPtr->_eRho<<" "
          //  <<setw(26)<<RhoeExt[iii]<<" "//exact phi, used for ehdplanner test, Lopez 2011 Table1
          <<setw(26)<<SPHPtPtr->_ePhi<<" "
          <<setw(26)<<SPHPtPtr->_eEx<<" "
          <<setw(26)<<SPHPtPtr->_eEy<<" "
          <<setw(26)<<SPHPtPtr->_NumNegbor<<" "//error e, used for ehdplanner test
          <<endl;
      }
    }

  //output time evolution of Q
  t=TimeSteps*Region._ControlSPH._DeltaT; 
 
  double R=0.05;
  double rhoe0=0.5;
  double Q=0.0;
  double Q0=0.0;
  CBasePt * PtPtr;
  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      if(Region._PtList[i]._Type==enSPHPt)
        {
          PtPtr=&Region._PtList[i];

          Q+=PtPtr->_mrho*PtPtr->_eRho;
          if(PtPtr->_PID==1)
            Q0+=PtPtr->_mrho*rhoe0;
        }
    
    }
  Q/=Q0;
  
  outfile.close();  
  outfilename.str("");
  outfilename<<"Q-evol"<<outputname<<".dat";
  outfile.open(outfilename.str().data(),ios::out|ios::app);
  if (!outfile)
    {
      cerr<<"Out put file 2  could not be open."<<endl;
      return;
    }
  outfile
    <<setiosflags(ios_base::scientific)
    <<setprecision(16)
    <<setw(26)<<t<<" "
    <<setw(26)<<Q<<" "
    <<endl;

  //output Er, from center to out
  outfile.close();
  double ErExact;
  outfilename.str("");
  outfilename<<"Er-"<<outputname<<TimeSteps<<".dat";
  outfile.open(outfilename.str().data(),ios::out);
  if (!outfile)
    {
      cerr<<"Out put file y=0 line could not be open."<<endl;
      return;
    }

  
  outfile
    <<"VARIABLES= "
    <<"\"AXISX\" "
    <<"\"Er\" "
    <<"\"ErExt\" "
    <<endl;
  outfile<<"ZONE T="<<"\""<<Region._ControlSPH._InfileName<<" t= "<<t<<"\" I= "<<34<<", F=POINT"<<endl;
  // for(unsigned int i=2244;i!=2278;++i) //N=64
  for(unsigned int i=8579;i!=8645;++i) //N=128
    {
      ErExact=(Region._PtList[i]._x<R?0:Q0/(2*PI*2*Region._PtList[i]._x));
      //if(ISZERO(Region._PtList[i]._y))
      {
        outfile
          <<setiosflags(ios_base::scientific)
          <<setprecision(16)
          <<setw(26)<<Region._PtList[i]._x<<" "
          <<setw(26)<<Region._PtList[i]._eEx<<" "
          <<setw(26)<<ErExact<<" "
          <<endl;
      }
    }

  

  RhoeExt.clear();

  cout<<"Output file "<<"\"" <<outfilename.str() <<"\""<<" has been written."<<endl;

  outfile.close();
}

void CKFile::outTecplotEHDPLANNER(CRegion & Region,unsigned int TimeSteps, string outputname)
{
  ofstream outfile;
  ostringstream outfilename;
  outfilename<<outputname<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".dat";
  outfile.open(outfilename.str().data(),ios::out);

  if (!outfile)
    {
      cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
      return;
    }


  outfile
    <<"VARIABLES= "
    <<"\"AXISX\" "
    <<"\"AXISY\" "
    //<<"\"AXISZ\" "
    <<"\"PID\" "
    <<"\"ID\" "
    <<"\"ID2\" "
    <<"\"DENSITY\" "
    //<<"\"DENSITY0\" "
    <<"\"VELOCITYX\" "
    <<"\"VELOCITYY\" "
    //<<"\"VELOCITYZ\" "
    <<"\"PRESSURE\" "
    <<"\"SMOOTHLENGTH\" "
    <<"\"PTSIZE\" "//particle size, for visualization simpilicity in tecplot
    <<"\"IFLAG\" "
    <<"\"DRHOE\" "
    <<"\"Fex\" "
    <<"\"Fey\" "
    <<"\"EPSILON\" "
    <<"\"KAPPA\" "
    <<"\"RHOE\" "
    <<"\"PHI\" "
    <<"\"EX\" "
    <<"\"EY\" "
    <<"\"PHIext\" "
    <<"\"Eext\" "
    <<"\"ErrorPHI\" "
    <<"\"ErrorE\" "
    <<"\"NUMNEIGHBOR\" "
      
    <<endl;
  outfile<<"ZONE T="<<"\""<<Region._ControlSPH._InfileName<<" step "<<TimeSteps<<"\" I= "<<Region._PtList.size()<<", F=POINT"<<endl;//I= 后面必须有个空格，否则重启读数不行
  CSPHPt *SPHPtPtr;

   //calculate the exact values of electric field and electric potential, and the errors between simulation and exact
  double beta,eta;//

  vector <double> Eext;
  vector <double> Phiext;
  vector <double> ErrorE;
  vector <double> ErrorPhi;
  Eext.resize(Region._PtList.size(),0.0);
  Phiext.resize(Region._PtList.size(),0.0);
  ErrorE.resize(Region._PtList.size(),0.0);
  ErrorPhi.resize(Region._PtList.size(),0.0);

  double y;
  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      SPHPtPtr=(CSPHPt *)(&Region._PtList[i]);

      if(SPHPtPtr->_Type==enSPHPt)
        {
          y=SPHPtPtr->_y;
          // //case 1: dielectric-dielectric
          // beta=3;
          // if(SPHPtPtr->_y<0)//lower
          //   {
          //     Phiext[i]=(-2*y+beta)/(1+beta);
          //     Eext[i]=2/(1+beta);
          //   }
          // else//upper
          //   {
          //     Phiext[i]=beta*(-2*y+1)/(1+beta);
          //     Eext[i]=2*beta/(1+beta);
          //   }

          //case 2: conducting-conducting
          beta=3;
          eta=2;
          //a=10;
          if(SPHPtPtr->_y<0)//lower
            {
              Phiext[i]=(-2*y+eta)/(1+eta);
              Eext[i]=2/(1+eta);
            }
          else//upper
            {
              Phiext[i]=eta*(-2*y+1)/(1+eta);
              Eext[i]=2*eta/(1+eta);
            }
          // //case 3: conducting-dielectric
          //      beta=3;
          // if(SPHPtPtr->_y<0)//lower
          //   {
          //     Phiext[i]=1;
          //     Eext[i]=0;
          //   }
          // else//upper
          //   {
          //     Phiext[i]=-2*y+1;
          //     Eext[i]=2;
          //   }

          //if case 3, lower, no errorE
          if(!ISZERO(Eext[i]))
            {
              ErrorE[i]=1-SPHPtPtr->_eEy/Eext[i];
            }
          ErrorPhi[i]=1-SPHPtPtr->_ePhi/Phiext[i];     
        }
    }
  

  for (size_t iii=0;iii<Region._PtList.size();iii++)
    {
      SPHPtPtr=(CSPHPt *)(&Region._PtList[iii]);
      {
        outfile
          <<setiosflags(ios_base::scientific)
          <<setprecision(16)
          <<setw(26)<<SPHPtPtr->_x<<" "
          <<setw(26)<<SPHPtPtr->_y<<" "
          //<<setw(26)<<SPHPtPtr->_z<<" "
          <<setw(26)<<SPHPtPtr->_PID<<" "
          <<setw(26)<<SPHPtPtr->_ID<<" "
          <<setw(26)<<SPHPtPtr->_ID2<<" "
          <<setw(26)<<SPHPtPtr->_rho<<" "
          //<<setw(26)<<SPHPtPtr->_rho0<<" "
          <<setw(26)<<SPHPtPtr->_u<<" "
          <<setw(26)<<SPHPtPtr->_v<<" "
          //<<setw(26)<<SPHPtPtr->_w<<" "
          <<setw(26)<<SPHPtPtr->_p<<" "
          <<setw(26)<<SPHPtPtr->_h<<" "
          <<setw(26)<<SPHPtPtr->_h/(Region._PartList[SPHPtPtr->_PID-1]._HdivDp)<<" "
          <<setw(26)<<SPHPtPtr->_Iflag<<" "
          <<setw(26)<<SPHPtPtr->_deRho<<" "
          <<setw(26)<<SPHPtPtr->_Fex<<" "
          <<setw(26)<<SPHPtPtr->_Fey<<" "
          <<setw(26)<<SPHPtPtr->_eEpsilon<<" "
          <<setw(26)<<SPHPtPtr->_eKappa<<" "
          <<setw(26)<<SPHPtPtr->_eRho<<" "
          <<setw(26)<<SPHPtPtr->_ePhi<<" "
          <<setw(26)<<SPHPtPtr->_eEx<<" "
          <<setw(26)<<SPHPtPtr->_eEy<<" "
          <<setw(26)<<Phiext[iii]<<" "//exact phi, used for ehdplanner test, Lopez 2011 Table1
          <<setw(26)<<Eext[iii]<<" "//exact e,used for ehdplanner test
          <<setw(26)<<ErrorPhi[iii]<<" "//error phi, used for ehdplanner test
          <<setw(26)<<ErrorE[iii]<<" "//error e, used for ehdplanner test
          <<setw(26)<<SPHPtPtr->_NumNegbor<<" "//error e, used for ehdplanner test
          <<endl;
      }
    }


  Phiext.clear();
  Eext.clear();
  ErrorE.clear();
  ErrorPhi.clear();
  
  cout<<"Output file "<<"\"" <<outfilename.str() <<"\""<<" has been written."<<endl;

  outfile.close();
}

//void CKFile::OutGrid(CExtGrid & ExtGrid,unsigned int TimeSteps)
//{
//	ofstream outfile;
//	ostringstream outfilename;
//	outfilename<<"GRID"<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".dat";
//	outfile.open(outfilename.str().data(),ios::out);
//	if (!outfile)
//	{
//		cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
//		return;
//	}
//			outfile
//			<<"VARIABLES= "
//			<<"\"AXISX\" "
//			<<"\"AXISY\" "
//			<<"\"AXISZ\" "
//			<<"\"CS\" "
//			<<"\"DIV\" "
//			<<endl;
//			outfile<<"ZONE T=\" step \" I="<<ExtGrid._GridNodeList.size()<<", F=POINT"<<endl;
//
//			CBaseGridNode *GridNode;
//			for (size_t iii=0;iii<ExtGrid._GridNodeList.size();iii++)
//	        {
//				GridNode=(CBaseGridNode *)ExtGrid._GridNodeList[iii];
//				outfile
//					<<setiosflags(ios_base::scientific)
//					<<setprecision(16)
//					<<setw(26)<<GridNode->_x<<" "
//					<<setw(26)<<GridNode->_y<<" "
//					<<setw(26)<<GridNode->_z<<" "
//					<<setw(26)<<GridNode->_C<<" "
//					<<setw(26)<<GridNode->_div<<" "
//					<<endl;
//	        }
//	outfile.close();
//}




void CKFile::Clear(CRegion &Region)
{
	Region._PtList.clear();
}

void CKFile::OutVector(std::vector<double> &v)//输出矩阵与向量，测试矩阵方程计算是否正确用
{
	ofstream outfile;
	ostringstream outfilename;
	outfilename<<"vector"<<".dat";
	outfile.open(outfilename.str().data(),ios::out);
	if (!outfile)
	{
		cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
		return;
	}

	for (size_t iii=0;iii<v.size();iii++)
  {
		outfile
			<<setiosflags(ios_base::scientific)
			<<setprecision(8)
			<<setw(16)<<v[iii]<<endl;
	}
	outfile.close();

}
void CKFile::OutMatrix(std::vector<CMatrix> &M,unsigned int dim)
{
	ofstream outfile;
	ostringstream outfilename;
	outfilename<<"matrix"<<".dat";
	outfile.open(outfilename.str().data(),ios::out);
	if (!outfile)
	{
		cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
		return;
	}

	double *a;
	a=(double *)malloc(dim*dim*sizeof(double));

	for(unsigned int i=0;i!=dim*dim;++i)
	{
		a[i]=0.0;
	}

	unsigned int rowid,colid;
	for(unsigned int i=0;i!=M.size();++i)
	{
		rowid=M[i]._RowID;
		colid=M[i]._ColID;

		a[(rowid-1)*dim+colid-1]=M[i]._Ele;
	}

	unsigned int id;
	for(unsigned int i=0;i!=dim;++i)
	{
		for(unsigned int j=0;j!=dim;++j)
		{
			id=i*dim+j;

			outfile
				<<setiosflags(ios_base::scientific)
				<<setprecision(5)
				<<setw(16)<<a[id]<<" ";
		}

		outfile<<endl;
	}
	outfile.close();
}
void CKFile::OutVTK(CRegion & Region,unsigned int TimeSteps, string outputname)
{
  ofstream outfile;
  ostringstream outfilename;
  outfilename<<outputname<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".vtk";
  outfile.open(outfilename.str().data(),ios::out);

  if (!outfile)
    {
      cerr<<"Out put file could not be open."<<endl<<endl;
      return;
    }

  unsigned int PtNum=Region._PtList.size();

  unsigned int icount=0;//used to create a new line
  unsigned int num1line=10;//numbers show in each line

  outfile<<"# vtk DataFile Version 2.0"<<endl;

  outfile<<outfilename.str()<<endl;

  outfile<<"ASCII"<<endl;

  outfile<<"DATASET UNSTRUCTURED_GRID"<<endl;

  outfile<<endl;

  outfile<<"POINTS "<<PtNum<<" double"<<endl;

  for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._x<<" "
         <<setw(26)<<Region._PtList[i]._y<<" "
         <<setw(26)<<0.0<<" "
         <<endl;
     }

   outfile<<"POINT_DATA "<<PtNum<<endl;

   outfile<<"SCALARS "<<"PID"<<" int"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile<<Region._PtList[i]._PID<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;
   
   outfile<<"SCALARS "<<"ID"<<" int"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile<<Region._PtList[i]._ID<<" ";

       icount++;
       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;
   outfile<<"SCALARS "<<"ID2"<<" int"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile<<Region._PtList[i]._ID2<<" ";
       
       icount++;
       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"IFLAG"<<" int"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile<<Region._PtList[i]._Iflag<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"NumberNeighbor"<<" int"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile<<Region._PtList[i]._NumNegbor<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;
   
   outfile<<"SCALARS "<<"RHO"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._rho<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;
   
   outfile<<"SCALARS "<<"PTSIZE"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._h/Region._PartList[Region._PtList[i]._PID-1]._HdivDp<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"H"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._h<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"Pressure"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._p<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"Colour"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._C<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"ElecPhi"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._ePhi<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"ElecEpsilon"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._eEpsilon<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"ElecKappa"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._eKappa<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"ElecRho"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._eRho<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   // <<setw(26)<<Phiext[iii]<<" "//exact phi, used for ehdplanner test, Lopez 2011 Table1
   // <<setw(26)<<Eext[iii]<<" "//exact e,used for ehdplanner test
   // <<setw(26)<<ErrorPhi[iii]<<" "//error phi, used for ehdplanner test
   // <<setw(26)<<ErrorE[iii]<<" "//error e, used for ehdplanner test

   //---------------------------------------------------------------------
   outfile<<"VECTORS "<<"Velocity"<<" double"<<endl;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._u
         <<setw(26)<<Region._PtList[i]._v
         <<setw(26)<<0.0
         <<endl;
     }
   
   outfile<<"VECTORS "<<"ElecForce"<<" double"<<endl;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._Fex
         <<setw(26)<<Region._PtList[i]._Fey
         <<setw(26)<<0.0
         <<endl;
     }

   outfile<<"VECTORS "<<"ElecE"<<" double"<<endl;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Region._PtList[i]._eEx
         <<setw(26)<<Region._PtList[i]._eEy
         <<setw(26)<<0.0
         <<endl;
     }

   //the following for calculating errors of ehd planner test
   double beta=Region._PartList[1]._eEpsilon/Region._PartList[2]._eEpsilon;//used for ehdplanner, lopez 2011 4.1 
   CBasePt * PtPtr;
   vector <double> Eext;
   vector <double> Phiext;
   vector <double> ErrorE;
   vector <double> ErrorPhi;
   Eext.resize(Region._PtList.size(),0.0);
   Phiext.resize(Region._PtList.size(),0.0);
   ErrorE.resize(Region._PtList.size(),0.0);
   ErrorPhi.resize(Region._PtList.size(),0.0);

   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       PtPtr=(&Region._PtList[i]);

       if(PtPtr->_Type==enSPHPt)
         {
           Eext[i]=(PtPtr->_PID-2)*2*beta/(1+beta)+(3-PtPtr->_PID)*2/(1+beta);
           Phiext[i]=(PtPtr->_PID-2)*beta*(-2*PtPtr->_y+1)/(1+beta)+(3-PtPtr->_PID)*(-2*PtPtr->_y+beta)/(1+beta);

           ErrorE[i]=1-PtPtr->_eEy/Eext[i];
           ErrorPhi[i]=1-PtPtr->_ePhi/Phiext[i];     
         }
     }
   
   outfile<<"SCALARS "<<"ElecEexact"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Eext[i]<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"ElecPhiexact"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<Phiext[i]<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;

   outfile<<"SCALARS "<<"ErrorElecE"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<ErrorE[i]<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;
   
   outfile<<"SCALARS "<<"ErrorElecPhi"<<" double"<<endl;
   outfile<<"LOOKUP_TABLE default"<<endl;
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       outfile
         <<setiosflags(ios_base::scientific)
         <<setprecision(16)
         <<setw(26)<<ErrorPhi[i]<<" ";
       icount++;

       if(icount%num1line==0)
         outfile<<endl;
     }
   outfile<<endl;
   
   Eext.clear();
   Phiext.clear();
   ErrorE.clear();
   ErrorPhi.clear();
   cout<<"vtk Output file "<<"\"" <<outfilename.str() <<"\""<<" has been written."<<endl;

   outfile.close();
}


void CKFile::InputRestart(CRegion & Region)//重启程序时的读模型部分
{
	string keyword;
	string line;
	vector<double> X,Y,Z;
	double value;
	unsigned int PID,SECID,EOSID,ID;
	unsigned int PtNum;
	char tch;

	ifstream infile;
	ostringstream infilename;
	infilename<<Region._ControlSPH._InfileName<<"TEC"<<setw(8)<<setfill('0')<<Region._ControlSPH._StartStep<<setw(4)<<".dat";
	infile.open(infilename.str().data(),ios::in);
	if (!infile)
	{
		cout<<"KFile could not be open."<<endl;
		return;
	}

	else
	{
		cout<<"Restart Step: "<<Region._ControlSPH._StartStep<<endl;
		cout<<"Start Reading Restart Model File."<<endl;
	}

	getline(infile,line);

	for(;;)
	{
		infile>>keyword;

		if(keyword=="I=")
		{
			infile>>value;
			PtNum=size_t(value);

			getline(infile,line);

			break;
		}
	}

	for(unsigned int i=0;i!=PtNum;i++)
	{
		CBasePt BasePt;

		infile>>value;
		BasePt._x=double(value);

		infile>>value;
		BasePt._y=double(value);

		infile>>value;
		BasePt._PID=size_t(value);

		infile>>value;
		BasePt._ID=size_t(value);

		infile>>value;
		BasePt._ID2=size_t(value);

		infile>>value;
		BasePt._rho=double(value);

		infile>>value;
		BasePt._u=double(value);

		infile>>value;
		BasePt._v=double(value);

		infile>>value;
		BasePt._p=double(value);

		infile>>value;
		BasePt._h=double(value);

		infile>>value;
		BasePt._Iflag=size_t(value);

		infile>>value;
		BasePt._VisEta=double(value);

		infile>>value;
		BasePt._VisGamma=double(value);

		infile>>value;
		BasePt._Nnx=double(value);

		infile>>value;
		BasePt._Nny=double(value);

		infile>>value;
		BasePt._AccBndx=double(value);

		infile>>value;
		BasePt._AccBndy=double(value);


		getline(infile,line);

		Region._PtList.push_back(BasePt);
	}

	infile.close();

	for(unsigned int i=0;i!=PtNum;++i)
	{
		PID=Region._PtList[i]._PID;

		Region._PtList[i]._r=Region._PtList[i]._h*2;

		Region._PtList[i]._Volume=pow(Region._PtList[i]._h/Region._PartList[PID-1]._HdivDp,2);

		Region._PtList[i]._rho0=_InitSPH[PID-1][2];
		Region._PtList[i]._T =_InitSPH[PID-1][3];

		Region._PtList[i]._m=Region._PtList[i]._rho0*Region._PtList[i]._Volume;

		Region._PtList[i]._C0=Region._PartList[PID-1]._C0;

		Region._PtList[i]._Cs=Region._ControlSPH._Cs;//声速赋初值

		if(Region._ControlSPH._SPHST==1)
		{
			Region._PtList[i]._H=Region._PartList[PID-1]._STHvsh*Region._PtList[i]._h;
			Region._PtList[i]._rr=2*Region._PtList[i]._H;
		}


		switch(Region._PartList[PID-1]._PartType)
		{
	 	  case enSPH:
					Region._PtList[i]._Type=enSPHPt;
					break;
	 	  case enNULL:
					Region._PtList[i]._Type=enNULLPt;
					break;
	 	  case enBnd:
					Region._PtList[i]._Type=enBndPt;
					break;
	 	  case enGhost1:
					Region._PtList[i]._Type=enGhost1Pt;
					break;
	 	  case enGhost2:
					Region._PtList[i]._Type=enGhost2Pt;
					break;
			case enDummy:
					Region._PtList[i]._Type=enDumPt;
					break;

			default:
					break;
		}

	}
}

void CKFile::InputMesh( CRegion & Region )
{
	string keyword;
	string line;
	vector<double> X,Y,Z;
	float value;
	unsigned int PID,SECID,EOSID,ID;
	char tch;

	ifstream infile;
	ostringstream infilename;
	infilename<<Region._ControlSPH._InfileName<<".msh";
	infile.open(infilename.str().data(),ios::in);
	if (!infile)
	{
		cout<<"Mesh File could not be open."<<endl;
		return;
	}

	else
	{
		cout<<"Start Reading Mesh File."<<endl;
	}


	while(!infile.eof())
	{
		infile>>keyword;

		if(keyword=="*MESH")
		{
			getline(infile,line);

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._MshMinX=double(value);

			infile>>value;
			infile.ignore();
			Region._ControlSPH._MshMaxX=double(value);

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._MshMinY=double(value);

			infile>>value;
			infile.ignore();
			Region._ControlSPH._MshMaxY=double(value);

			infile>>value;
			//infile.ignore();
			Region._ControlSPH._MshSize=double(value);

		}

		if(keyword=="*NODE")
		{
			for(;;)
			{
				infile>>tch;
				infile.seekg(-1,ios::cur);

				if((tch=='*')||(tch=='$'))
				{
					break;
				}

				CNode TempNode;

				infile>>value;
				infile.ignore();
				TempNode._ID=size_t(value);

				infile>>value;
				infile.ignore();
				TempNode._x=double(value);

				infile>>value;
				infile.ignore();
				TempNode._y=double(value);

				for(int i=0;i<3;i++)
				{
					infile>>value;
					infile.ignore();
				}

				Region._NodeList.push_back(TempNode);
			}
		}

		if(keyword=="*ELEMENT_SHELL_THICKNESS")
		{
			for(;;)
			{
				vector<double>tPosx(4);
				vector<double>tPosy(4);
				vector<double>tPartile(3);

				CMesh TempMesh;

				infile>>tch;
				infile.seekg(-1,ios::cur);
				if((tch=='*')||(tch=='$'))	break;

				infile>>value;
				infile.ignore();
				ID=size_t(value);
				TempMesh._ID=ID;

				infile>>value;
				infile.ignore();
				PID=size_t(value);
				TempMesh._PID=PID;


				infile>>value;
				infile.ignore();
				ID=size_t(value);
				TempMesh._MeshNodList[0]=&Region._NodeList[ID-1];
				tPosx[0]=Region._NodeList[ID-1]._x;
				tPosy[0]=Region._NodeList[ID-1]._y;

				infile>>value;
				infile.ignore();
				ID=size_t(value);
				TempMesh._MeshNodList[1]=&Region._NodeList[ID-1];
				tPosx[1]=Region._NodeList[ID-1]._x;
				tPosy[1]=Region._NodeList[ID-1]._y;

				infile>>value;
				infile.ignore();
				ID=size_t(value);
				TempMesh._MeshNodList[2]=&Region._NodeList[ID-1];
				tPosx[2]=Region._NodeList[ID-1]._x;
				tPosy[2]=Region._NodeList[ID-1]._y;

				infile>>value;
				infile.ignore();
				ID=size_t(value);
				TempMesh._MeshNodList[3]=&Region._NodeList[ID-1];
				tPosx[3]=Region._NodeList[ID-1]._x;
				tPosy[3]=Region._NodeList[ID-1]._y;

				for(int i=0;i<4;i++)
				{
					infile>>value;
					infile.ignore();
				}

				Center(tPosx,tPosy,tPartile);

				TempMesh._x=tPartile[0];
				TempMesh._y=tPartile[1];
				TempMesh._Volume=tPartile[2];


				TempMesh._h=/*2.0**/Region._PartList[PID-1]._HdivDp*pow(Region._PtList[0]._Volume,0.5);//折衷，有点问题,将网格的光滑长度定义为粒子的光滑长度的两倍
				TempMesh._r=2.0*TempMesh._h;
				TempMesh._C0=1;
				TempMesh._Iflag=0;//对于背景网格，_Iflag代表其是否与SPH粒子发生作用

				TempMesh._H=Region._PartList[PID-1]._STHvsh*TempMesh._h;
				TempMesh._rr=2*TempMesh._H;

				//BasePt._Cs=Region._ControlSPH._Cs;//为每个粒子的声速赋初值
				//BasePt._u=_InitSPH[PID-1][0];
				//BasePt._v=_InitSPH[PID-1][1];

				//BasePt._rho0=Region._PtList[0]._rho0;//折衷，有问题，如果两个Part密度不同，则会出问题
				//BasePt._rho =_InitSPH[PID-1][2];
				//BasePt._T =_InitSPH[PID-1][3];

				//BasePt._m=BasePt._rho0*BasePt._Volume;



				//TempMesh._Type=enSPHPt;

				//BasePt._p=0;

				Region._MeshList.push_back(TempMesh);

				tPosx.clear();
				tPosy.clear();
				tPartile.clear();
			}
		}

		if(keyword=="*END")
		{
			break;
		}
	}

}

void CKFile::OutTecplotMsh( CRegion & Region,unsigned int TimeSteps )
{
	ofstream outfile;
	ostringstream outfilename;	
	outfilename<<Region._ControlSPH._InfileName<<"_PTMESH_"<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".dat";	
	outfile.open(outfilename.str().data(),ios::out);

	if (!outfile)
	{
		cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
		return;
	}
	outfile
		<<"VARIABLES= "
		<<"\"AXISX\" "
		<<"\"AXISY\" "
		//<<"\"AXISZ\" "
		<<"\"PID\" "
		<<"\"ID\" "
		<<"\"ID2\" "
		<<"\"DENSITY\" "
		//<<"\"DENSITY0\" "
		<<"\"VELOCITYX\" "
		<<"\"VELOCITYY\" "
		//<<"\"VELOCITYZ\" "
		<<"\"PRESSURE\" "
		<<"\"SMOOTHLENGTH\" "
		<<"\"IFLAG\" "
		<<"\"VISETA\" "
		<<"\"VISGAMMA\" "
		<<"\"DIV\" "
		<<"\"NX\" "
		<<"\"NY\" "
		<<"\"NnX\" "
		<<"\"NnY\" "
		<<"\"LAMBDA\" "
		<<"\"C0\" "
		<<"\"C\" "
		<<endl;
	outfile<<"ZONE T="<<"\" Mesh step "<<TimeSteps<<"\" I= "<<Region._MeshList.size()<<", F=POINT"<<endl;//I= 后面必须有个空格，否则重启读数不行

	CMesh * MshPtr;
	for (size_t iii=0;iii<Region._MeshList.size();iii++)
	{
		MshPtr=&Region._MeshList[iii];
		{	
			outfile
				<<setiosflags(ios_base::scientific)
				<<setprecision(16)
				<<setw(26)<<MshPtr->_x<<" "
				<<setw(26)<<MshPtr->_y<<" "
				//<<setw(26)<<MshPtr->_z<<" "
				<<setw(26)<<MshPtr->_PID<<" "
				<<setw(26)<<MshPtr->_ID<<" "
				<<setw(26)<<MshPtr->_ID2<<" "
				<<setw(26)<<0.0<<" "
				//<<setw(26)<<MshPtr->_rho0<<" "
				<<setw(26)<<0.0<<" "
				<<setw(26)<<0.0<<" "
				//<<setw(26)<<MshPtr->_w<<" "
				<<setw(26)<<0.0<<" "
				<<setw(26)<<MshPtr->_h<<" "
				<<setw(26)<<MshPtr->_Iflag<<" "
				<<setw(26)<<0.0<<" "
				<<setw(26)<<0.0<<" "
				<<setw(26)<<MshPtr->_Curvature<<" "
				<<setw(26)<<MshPtr->_nx<<" "
				<<setw(26)<<MshPtr->_ny<<" "
				<<setw(26)<<MshPtr->_Nnx<<" "
				<<setw(26)<<MshPtr->_Nny<<" "
				<<setw(26)<<0.0<<" "
				<<setw(26)<<MshPtr->_C0<<" "
				<<setw(26)<<MshPtr->_C<<" "
				<<endl;
		}
	}


	cout<<"Output mesh file has been written."<<endl;

	outfile.close();

}
void CKFile::OutMsh( CRegion & Region,unsigned int TimeSteps,double CurTime )//按网格形式输出网格信息
{
	ofstream outfile;
	ostringstream outfilename;	
	outfilename<<Region._ControlSPH._InfileName<<"_MESH_"<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".dat";	
	outfile.open(outfilename.str().data(),ios::out);

	if (!outfile)
	{
		cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
		return;
	}

	outfile<<"TITLE="<<"\"SPH BASE MESH\""<<endl;

	outfile<<"FILETYPE=FULL"<<endl;

	outfile<<"VARIABLES= "
			 	<<"\"AXISX\" "
				<<"\"AXISY\" "
				//<<"\"AXISZ\" "
				<<"\"PID\" "
				<<"\"ID\" "
				<<"\"ID2\" "
				<<"\"DENSITY\" "
				//<<"\"DENSITY0\" "
				<<"\"VELOCITYX\" "
				<<"\"VELOCITYY\" "
				//<<"\"VELOCITYZ\" "
				<<"\"PRESSURE\" "
				<<"\"SMOOTHLENGTH\" "
				<<"\"IFLAG\" "
				<<"\"VISETA\" "
				<<"\"VISGAMMA\" "
				<<"\"DIV\" "
				<<"\"NX\" "
				<<"\"NY\" "
				<<"\"NnX\" "
				<<"\"NnY\" "
				<<"\"LAMBDA\" "
				<<"\"C0\" "
				<<"\"C\" "
				<<endl;

	outfile<<"ZONE"<<endl;

	outfile<<"T="<<"\""<<Region._ControlSPH._InfileName<<"_MESH_"<<TimeSteps<<"\" "<<endl;

	outfile<<"STRANDID="<<TimeSteps<<", "<<"SOLUTIONTIME="<<CurTime<<endl;

	outfile<<"ZONETYPE=FEQUADRILATERAL"<<endl;

	outfile<<"NODES="<<Region._NodeList.size()<<","<<"ELEMENTS="<<Region._MeshList.size()<<endl;

	outfile<<"FACENEIGHBORMODE=LOCALONETOONE"<<endl;

	outfile<<"DATAPACKING=BLOCK"<<endl;

	outfile<<"VARLOCATION=([3-21]=CELLCENTERED)"<<endl;
	
	outfile<<setiosflags(ios_base::scientific)<<setprecision(16);
		
	//1,2 output nodal variables x,y
	for (unsigned int i=0;i!=Region._NodeList.size();++i)
	{
		outfile/*<<setiosflags(ios_base::scientific)
					<<setprecision(16)*/
					<<setw(26)<<Region._NodeList[i]._x<<" ";
		
		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;

	for (unsigned int i=0;i!=Region._NodeList.size();++i)
	{
		outfile/*<<setiosflags(ios_base::scientific)
					<<setprecision(16)*/
					<<setw(26)<<Region._NodeList[i]._y<<" ";
		
		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;

	//output cell-centered variables
	//3.PID
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<1<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;	
	
	//4.ID
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._ID<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	
	//5.ID2
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._ID2<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;

	//6.rho
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<0.0<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	//7.u
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<0.0<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	//8.v
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<0.0<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	//9.p
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<0.0<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;

	//10.h
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._h<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;

	//11.Iflag
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._Iflag<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	//12.eta
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<0.0<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	//13.gamma
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<0.0<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	
	//14.curvature
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._Curvature<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;

	//15.normal x
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._nx<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	//16.normal y
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._ny<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	//17.unit normal x
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._Nnx<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;
	//18.unit normal y
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._Nny<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;

	//19.lambda
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<0.0<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;

	//20.colour
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._C0<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;

	//21.interpolated colour
	for (unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile
			<<setw(26)<<Region._MeshList[i]._C<<" ";

		if((i+1)%20==0)
			outfile<<endl;
	}
	outfile<<endl;


	//output element nodes numbers
	for(unsigned int i=0;i!=Region._MeshList.size();++i)
	{
		outfile<<Region._MeshList[i]._MeshNodList[0]->_ID<<" "
			    <<Region._MeshList[i]._MeshNodList[1]->_ID<<" "
	        <<Region._MeshList[i]._MeshNodList[2]->_ID<<" "
	        <<Region._MeshList[i]._MeshNodList[3]->_ID<<" "
					<<endl;
	}

	outfile.close();
}

//for test, output particle pairs,2019.11.01
void CKFile::OutPtNeighbour(CRegion & Region, CBasePt * PtPtr,unsigned int TimeSteps, string outputname)
{
  CBasePt * PtiPtr,*PtjPtr;
  CKnl * KnlPtr;
  unsigned int ID;

  vector<CBasePt *> NeighborList;
  NeighborList.clear();

  //the first the neighbor if ptptr itself
  NeighborList.push_back(PtPtr);

  //find the neighbours of ptptr
  for(unsigned int i=0;i!=Region._PtPairList.size();++i)
    {
      PtiPtr=Region._PtPairList[i]._PtiPtr;
      PtjPtr=Region._PtPairList[i]._PtjPtr;

      if(PtiPtr!=PtjPtr)
        {          
          if(PtiPtr==PtPtr)
            NeighborList.push_back(PtjPtr);

          if(PtjPtr==PtPtr)
            NeighborList.push_back(PtiPtr);
        }
    }

  ofstream outfile;
  ostringstream outfilename;
  outfilename<<outputname<<setw(8)<<setfill('0')<<TimeSteps<<setw(4)<<".dat";
  outfile.open(outfilename.str().data(),ios::out);

  if (!outfile)
    {
      cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
      return;
    }

  outfile
    <<"VARIABLES= "
    <<"\"AXISX\" "
    <<"\"AXISY\" "
    //<<"\"AXISZ\" "
    <<"\"PID\" "
    <<"\"ID\" "
    <<"\"ID2\" "
    <<"\"DENSITY\" "
    //<<"\"DENSITY0\" "
    <<"\"VELOCITYX\" "
    <<"\"VELOCITYY\" "
    //<<"\"VELOCITYZ\" "
    <<"\"PRESSURE\" "
    <<"\"SMOOTHLENGTH\" "
    <<"\"PTSIZE\" "//particle size, for visualization simpilicity in tecplot
    <<"\"IFLAG\" "
    <<"\"CI\" "
    <<"\"Fex\" "
    <<"\"Fey\" "
    <<"\"EPSILON\" "
    <<"\"KAPPA\" "
    <<"\"RHOE\" "
    <<"\"PHI\" "
    <<"\"EX\" "
    <<"\"EY\" "
    // <<"\"PHIext\" "
    // <<"\"Eext\" "
    // <<"\"ErrorPHI\" "
    // <<"\"ErrorE\" "
     <<"\"NUMNEIGHBOR\" "
      
    <<endl;
  outfile<<"ZONE T="<<"\""<<"NeighboursofPt"<<PtPtr->_ID<<" step "<<TimeSteps<<"\" I= "<<NeighborList.size()<<", F=POINT"<<endl;//I= 后面必须有个空格，否则重启读数不行
  CSPHPt *SPHPtPtr;

  for (size_t iii=0;iii!=NeighborList.size();++iii)
    {
      SPHPtPtr=(CSPHPt *)(NeighborList[iii]);
      {
        outfile
          <<setiosflags(ios_base::scientific)
          <<setprecision(16)
          <<setw(26)<<SPHPtPtr->_x<<" "
          <<setw(26)<<SPHPtPtr->_y<<" "
          //<<setw(26)<<SPHPtPtr->_z<<" "
          <<setw(26)<<SPHPtPtr->_PID<<" "
          <<setw(26)<<SPHPtPtr->_ID<<" "
          <<setw(26)<<SPHPtPtr->_ID2<<" "
          <<setw(26)<<SPHPtPtr->_rho<<" "
          //<<setw(26)<<SPHPtPtr->_rho0<<" "
          <<setw(26)<<SPHPtPtr->_u<<" "
          <<setw(26)<<SPHPtPtr->_v<<" "
          //<<setw(26)<<SPHPtPtr->_w<<" "
          <<setw(26)<<SPHPtPtr->_p<<" "
          <<setw(26)<<SPHPtPtr->_h<<" "
          <<setw(26)<<SPHPtPtr->_h/(Region._PartList[SPHPtPtr->_PID-1]._HdivDp)<<" "
          <<setw(26)<<SPHPtPtr->_Iflag<<" "
          <<setw(26)<<SPHPtPtr->_C<<" "
          <<setw(26)<<SPHPtPtr->_Fex<<" "
          <<setw(26)<<SPHPtPtr->_Fey<<" "
          <<setw(26)<<SPHPtPtr->_eEpsilon<<" "
          <<setw(26)<<SPHPtPtr->_eKappa<<" "
          <<setw(26)<<SPHPtPtr->_eRho<<" "
          <<setw(26)<<SPHPtPtr->_ePhi<<" "
          <<setw(26)<<SPHPtPtr->_eEx<<" "
          <<setw(26)<<SPHPtPtr->_eEy<<" "
          // <<setw(26)<<Phiext[iii]<<" "//exact phi, used for ehdplanner test, Lopez 2011 Table1
          // <<setw(26)<<Eext[iii]<<" "//exact e,used for ehdplanner test
          // <<setw(26)<<ErrorPhi[iii]<<" "//error phi, used for ehdplanner test
          // <<setw(26)<<ErrorE[iii]<<" "//error e, used for ehdplanner test
         <<setw(26)<<SPHPtPtr->_NumNegbor<<" "//error e, used for ehdplanner test
          <<endl;
      }
    }


  cout<<"Details of neighbours of particle "<<PtPtr->_ID
      <<"has been written in file "<<"\"" <<outfilename.str() <<"\""<<endl;

  NeighborList.clear();
  outfile.close();
  
}
