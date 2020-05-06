#include "CalBndNorm.h"

CCalBndNorm::CCalBndNorm()
{
}

CCalBndNorm::~CCalBndNorm()
{

}

void CCalBndNorm::Solve(CRegion & Region)
{
	CBasePt * PtPtr;

	RegionBnd._ControlSPH._InfileName=Region._ControlSPH._InfileName;
	
	InputBndAssit(RegionBnd);

	GetNblBnd(RegionBnd);

	GetNormBnd(RegionBnd);

	//����������ı߽�ķ���ֵ��Region�еı߽�����
	for (unsigned int i=0;i!=Region._PtList.size();++i)
	{
			if (Region._PtList[i]._Type==enDumPt/*&&Region._PtList[i]._PID==1*/)
			{
				PtPtr=&Region._PtList[i];

				PtPtr->_Nnx=RegionBnd._PtList[i]._Nnx;
				PtPtr->_Nny=RegionBnd._PtList[i]._Nny;
			}

			//if (Region._PtList[i]._Type==enDummy&&Region._PtList[i]._PID==2)//control disk����
			//{
			//	PtPtr=&Region._PtList[i];

			//	if (PtPtr->_y>0.001)
			//	{
			//		PtPtr->_Nnx=0;
			//		PtPtr->_Nny=1;
			//	}

			//	else
			//	{
			//		PtPtr->_Nnx=0;
			//		PtPtr->_Nny=-1;
			//	}
			//}


		//����Ϊ����ģ�͵ķ���ֵ
	//	if (Region._PtList[i]._Type==enDummy&&(Region._PtList[i]._PID==1||Region._PtList[i]._PID==4))
	//	{
	//		PtPtr=&Region._PtList[i];

	//		PtPtr->_Nnx=RegionBnd._PtList[i]._Nnx;
	//		PtPtr->_Nny=RegionBnd._PtList[i]._Nny;

	//		if (PtPtr->_PID==1&&PtPtr->_x<0.005)
	//		{
	//			PtPtr->_Nnx=0.0;
	//			PtPtr->_Nny=1.0;
	//		}
	//	}

	//	if (Region._PtList[i]._Type==enDummy&&Region._PtList[i]._PID==2)//����ģ��
	//	{
	//		PtPtr=&Region._PtList[i];

	//		PtPtr->_Nnx=-RegionBnd._PtList[i]._Nnx;
	//		PtPtr->_Nny=-RegionBnd._PtList[i]._Nny;
	//	}
	//	if (Region._PtList[i]._Type==enDummy&&Region._PtList[i]._PID==3)//����ģ��
	//	{
	//		PtPtr=&Region._PtList[i];

	//		PtPtr->_Nnx=-RegionBnd._PtList[i]._Nnx;
	//		PtPtr->_Nny=-RegionBnd._PtList[i]._Nny;
	//	}

	//	if (Region._PtList[i]._Type==enDummy&&Region._PtList[i]._PID==4)//mould���������б��
	//	{
	//		PtPtr=&Region._PtList[i];

	//		if(PtPtr->_y>0.0)
	//		{
	//			PtPtr->_Nnx=0.9985;
	//			PtPtr->_Nny=0.0545;
	//		}
	//	}

	//	if (Region._PtList[i]._Type==enDummy&&(Region._PtList[i]._PID==1||Region._PtList[i]._PID==2))
	//	{
	//		PtPtr=&Region._PtList[i];

	//		if (PtPtr->_x<0.046&&PtPtr->_x>0.04
	//			&&PtPtr->_y<0.08&&PtPtr->_y>0.073)
	//		{
	//			PtPtr->_Nnx=1.0;
	//			PtPtr->_Nny=0.0;
	//		}

	//	}
	//---------------------------------------------------------------


		
		//if (Region._PtList[i]._Type==enDummy&&Region._PtList[i]._PID==2)//�����������е��˶����ڲ����ӵķ���
		//{
		//	PtPtr=&Region._PtList[i];

		//	PtPtr->_Nnx=-1*RegionBnd._PtList[i]._Nnx;
		//	PtPtr->_Nny=-1*RegionBnd._PtList[i]._Nny;
		//}

		//if (Region._PtList[i]._Type==enDumPt&&Region._PtList[i]._PID==2)//part 2,����
		//{
		//	PtPtr=&Region._PtList[i];

		//	PtPtr->_Nnx=1.0;
		//	PtPtr->_Nny=0.0;
		//}
	}

	RegionBnd._PtList.clear();
	RegionBnd._PtPairList.clear();
	RegionBnd._KnlList.clear();

}


void CCalBndNorm::InputBndAssit( CRegion & Region )
{
	string keyword;
	string line;
	vector<double> X,Y,Z;
	float value;
	unsigned int PID,SECID,EOSID,ID;
	char tch;

	ifstream infile;
	ostringstream infilename;
	infilename<<Region._ControlSPH._InfileName<<".bnd";	
	infile.open(infilename.str().data(),ios::in);
	if (!infile)
	{
		cout<<"KFile could not be open."<<endl;
		return;
	}

	else
	{
		cout<<"Start Reading Boundary Normal Assistant File."<<endl;
	}

	while(!infile.eof())
	{
		infile>>keyword;

		if(keyword=="*NODE")
		{
			if(Region._ControlSPH._StartStep==0)//��0ʱ�䲽��ʼ����ʱ����ò�������
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
			if(Region._ControlSPH._StartStep==0)//��0ʱ�䲽��ʼ����ʱ����ò�������
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

					CenterBnd(tPosx,tPosy,tPartile);

					BasePt._x=tPartile[0];
					BasePt._y=tPartile[1];
					BasePt._Volume=tPartile[2];

					BasePt._rho0=1000;
					BasePt._rho =1000;

					BasePt._m=BasePt._rho0*BasePt._Volume;

					BasePt._h=2*pow(BasePt._Volume,0.5);
					BasePt._r=2.0*BasePt._h;
					BasePt._C0=1;
					BasePt._C=0;

					BasePt._Cs=Region._ControlSPH._Cs;//Ϊÿ�����ӵ����ٸ���ֵ

					if(Region._ControlSPH._SPHST==1)
					{
						BasePt._H=Region._PartList[PID-1]._STHvsh*BasePt._h;
						BasePt._rr=2*BasePt._H;
					}

					BasePt._p=0;
					BasePt._Iflag=1;
					BasePt._Type=enSPHPt;								
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

}

void CCalBndNorm::CenterBnd( vector<double> &px,vector<double> &py,vector<double> &p )
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
	p[0]=(xc1*A1+xc2*A2)/A;//�����������x	
	p[1]=(yc1*A1+yc2*A2)/A;//�����������y
	p[2]=A;//������	

}

void CCalBndNorm::GetNblBnd( CRegion & Region )
{
	double dxiac;
	double dyiac;
	double dxiac2;
	double dyiac2;
	double driac2;
	double r,rr;
	double r2,rr2;
	unsigned int ID;
	CPtPair PtPair;

	//set min & max point
	double dp;

	//�ҵ���һ�����ڼ������е����ӣ��������ӻ�Null���ӣ�
	for(unsigned int i=0;i!=Region._PtList.size();++i)
	{
		if(Region._PtList[i]._Iflag==1)
		{
			_MaxPointx=Region._PtList[i]._x;
			_MaxPointy=Region._PtList[i]._y;

			_MinPointx=Region._PtList[i]._x;
			_MinPointy=Region._PtList[i]._y;

			_MinPointr=Region._PtList[i]._r;

			break;
		}
	}

	//�ҵ��������ڵ���С�����ֵ
	for(size_t i=0;i!=Region._PtList.size();++i)
	{
		if(Region._PtList[i]._Iflag==1)
		{
			if(_MaxPointx<Region._PtList[i]._x)	_MaxPointx=Region._PtList[i]._x;
			if(_MaxPointy<Region._PtList[i]._y)	_MaxPointy=Region._PtList[i]._y;

			if(_MinPointx>Region._PtList[i]._x)	_MinPointx=Region._PtList[i]._x;
			if(_MinPointy>Region._PtList[i]._y)	_MinPointy=Region._PtList[i]._y;

			if(_MinPointr>Region._PtList[i]._r)	_MinPointr=Region._PtList[i]._r;
		}
	}

	_MinPointr=Region._PtList[0]._r;//��ʱ�ģ��϶�������

	_MaxPointx+=_MinPointr*0.5;
	_MinPointx-=_MinPointr*0.5;

	_MaxPointy+=_MinPointr*0.5;
	_MinPointy-=_MinPointr*0.5;

	//��̬����ÿ�����ĺ�����
	dp=pow(Region._PtList[0]._Volume,0.5);//��ʱ

	_CellNumx=size_t((_MaxPointx-_MinPointx)/dp*1.0);
	_CellNumy=size_t((_MaxPointy-_MinPointy)/dp*1.0);


	cout<<"cell num x:"<<_CellNumx<<endl;
	cout<<"cell num y:"<<_CellNumy<<endl;

	_CellSizex=(_MaxPointx-_MinPointx)/double(_CellNumx);
	_CellSizey=(_MaxPointy-_MinPointy)/double(_CellNumy);

	//���ú����������ڲ����
	_CellList.resize(_CellNumx*_CellNumy,std::vector<CBasePt *>(0,NULL));

	for(size_t i=0;i<_CellList.size();i++)
	{
		_CellList[i].clear();
	}

	//PutPt2Cell(CRegion & Region)
	{
		double xM,yM;	
		size_t xI,yI;
		unsigned int ID;


		for(size_t i=0; i!=Region._PtList.size();i++)
		{
			if(Region._PtList[i]._Iflag==1)
			{
				xM=Region._PtList[i]._x-_MinPointx;
				yM=Region._PtList[i]._y-_MinPointy;

				xI=size_t(floor(xM/_CellSizex));
				yI=size_t(floor(yM/_CellSizey));

				_CellList[yI*_CellNumx+xI].push_back(&Region._PtList[i]);	
			}
		}	
	}

	//���������ӶԲ���
	//ֻ������������Ϊ��������
	for(size_t f=0;f<Region._PtList.size();f++)
	{
		ID=f;

		if(Region._PtList[ID]._Type==enSPHPt)
		{
			PrimarySchBnd(Region._PtList[ID]._x,Region._PtList[ID]._y,Region._PtList[ID]._r,_PtNbList);

			for(size_t i=0;i<_PtNbList.size();i++)//for all the number of the _PtNbList
			{
				dxiac=_PtNbList[i]->_x-Region._PtList[ID]._x;
				dyiac=_PtNbList[i]->_y-Region._PtList[ID]._y;

				dxiac2=dxiac*dxiac;
				dyiac2=dyiac*dyiac;

				driac2=dxiac2+dyiac2;
				r=_PtNbList[i]->_r+Region._PtList[ID]._r;
				r2=r*r;

				//the second searching
				if(4*driac2<=r2)
				{						
					switch (_PtNbList[i]->_Type)
					{
					case enSPHPt: 
						{	
							if((ID+1)<=_PtNbList[i]->_ID)
							{
								PtPair._Type=enSPHPtPair;
								PtPair._PtiPtr=&Region._PtList[ID];
								PtPair._PtjPtr=_PtNbList[i];
								PtPair._driac2=driac2;

								Region._PtPairList.push_back(PtPair);

								break;
							}

							else
							{
								double PtNbr;
								PtNbr=_PtNbList[i]->_r*_PtNbList[i]->_r;
								if(driac2>PtNbr)
								{
									PtPair._Type=enSPHPtPair;
									PtPair._PtiPtr=&Region._PtList[ID];
									PtPair._PtjPtr=_PtNbList[i];
									PtPair._driac2=driac2;

									Region._PtPairList.push_back(PtPair);
								}
								break;
							}
						}

					default:
						break;
					}

				}
			}
		}
	}


	//����˺���
	double icount;
	double x,y,h,distance2;
	CBasePt * PtiPtr,*PtjPtr;
	Region._KnlList.resize(Region._PtPairList.size());

	icount=0;
	for(size_t i=0;i<Region._PtPairList.size();i++)
	{
		{
			PtiPtr=Region._PtPairList[i]._PtiPtr;
			PtjPtr=Region._PtPairList[i]._PtjPtr;

			x=PtiPtr->_x-PtjPtr->_x;
			y=PtiPtr->_y-PtjPtr->_y;
			h=0.5*(PtiPtr->_h+PtjPtr->_h); 


			distance2=Region._PtPairList[i]._driac2;

			GetWBnd(x,y,h,distance2,Region._KnlList[icount]._W,Region._KnlList[icount]._Wx,Region._KnlList[icount]._Wy,Region._KnlList[icount]._Ww);

			icount++;
		}		
	}

	_CellList.clear();
	
}

void CCalBndNorm::PrimarySchBnd( double x, double y, double r,std::vector<CBasePt *> & PtNbList )
{
	double dxiac;
	double dyiac;
	double driac2;

	double dxiac2;
	double dyiac2;
	double r2;
	unsigned int ID;
	r2=r*r;

	PtNbList.resize(0,NULL);

	//determate the given pt to which cell
	double xM,yM;
	size_t xI,yI;
	xM=x-_MinPointx;
	yM=y-_MinPointy;

	xI=size_t(floor(xM/_CellSizex));
	yI=size_t(floor(yM/_CellSizey));

	//determate searching range
	int NxI;
	int NyI;
	int rsizex;
	int rsizey;

	rsizex=size_t(floor(r/_CellSizex)+1);
	rsizey=size_t(floor(r/_CellSizey)+1);

	for(int i=-rsizex;i<=rsizex;i++)
	{
		NxI=xI+i;
		if(NxI>=0&&NxI<_CellNumx)
		{
			for(int j=-rsizey;j<=rsizey;j++)
			{
				NyI=yI+j;
				if(NyI>=0&&NyI<_CellNumy)
				{
					for(size_t f=0;f<_CellList[NyI*_CellNumx+NxI].size();f++)
					{
						dxiac=_CellList[NyI*_CellNumx+NxI][f]->_x-x;
						dyiac=_CellList[NyI*_CellNumx+NxI][f]->_y-y;

						dxiac2=dxiac*dxiac;
						dyiac2=dyiac*dyiac;

						driac2=dxiac2+dyiac2;

						if(driac2<=r2)
						{ 
							PtNbList.push_back(_CellList[NyI*_CellNumx+NxI][f]);
						}
					}
				}
			}
		}
	}
}

void CCalBndNorm::GetNormBnd( CRegion & Region )
{
	CBasePt * PtiPtr,* PtjPtr;
	CBasePt * BasePtr;
	CKnl * KnlPtr;

	double mi,mj;
	double rhoi,rhoj;
	double Nnxx,Nnyy;
	unsigned int C0i,C0j;
	unsigned int IDi,IDj;
	unsigned int ID,ID2;
	unsigned int PID;
	double Axx,Ayx,Axy,Ayy;
	double xij,yij;
	double normx,normy;

	vector<double> bx;//CSPM������Դ��
	vector<double> by;
	vector<unsigned int> N;//�����ж��Ƿ�Լ��㷨����Ч

	bx.resize(Region._PtList.size(),0.0);
	by.resize(Region._PtList.size(),0.0);
	N.resize(Region._PtList.size(),0);

	for (unsigned int i=0;i!=Region._PtList.size();++i)
	{
		BasePtr=&Region._PtList[i];

		BasePtr->_CSPMAxx=0.0;
		BasePtr->_CSPMAxy=0.0;
		BasePtr->_CSPMAyx=0.0;
		BasePtr->_CSPMAyy=0.0;
	}

	//1.��ֵ�õ�ɫֵ
	for(unsigned int i=0;i!=Region._PtPairList.size();++i)
	{
		PtiPtr=Region._PtPairList[i]._PtiPtr;
		PtjPtr=Region._PtPairList[i]._PtjPtr;
		KnlPtr=&Region._KnlList[i];

		C0i=1;
		C0j=1;

		PtiPtr->_C+=PtjPtr->_m/PtjPtr->_rho*KnlPtr->_W*C0j;

		if(PtiPtr!=PtjPtr)
		{
			PtjPtr->_C+=PtiPtr->_m/PtiPtr->_rho*KnlPtr->_W*C0i;
		}
	}

	//2 ���㷨�������ķ���
	//2.1����CSPM����ϵ��
	for (size_t i=0;i<Region._PtPairList.size();i++)
	{
		//if(Region._PtPairList[i]._Type==enSPHPtPair||Region._PtPairList[i]._Type==enSPHNULLPtPair)
		{
			PtiPtr=Region._PtPairList[i]._PtiPtr;
			PtjPtr=Region._PtPairList[i]._PtjPtr;
			KnlPtr=&Region._KnlList[i];

			mi=PtiPtr->_m;
			mj=PtjPtr->_m;

			rhoi=PtiPtr->_rho;
			rhoj=PtjPtr->_rho;

			xij=PtiPtr->_x-PtjPtr->_x;
			yij=PtiPtr->_y-PtjPtr->_y;

			normx=mj/rhoj*KnlPtr->_Wx;
			normy=mj/rhoj*KnlPtr->_Wy;

			PtiPtr->_CSPMAxx+=xij*normx;
			PtiPtr->_CSPMAyx+=yij*normx;

			PtiPtr->_CSPMAxy+=xij*normy;
			PtiPtr->_CSPMAyy+=yij*normy;

			if(PtiPtr!=PtjPtr)
			{			
				normx=mi/rhoi*KnlPtr->_Wx;
				normy=mi/rhoi*KnlPtr->_Wy;

				PtjPtr->_CSPMAxx+=xij*normx;
				PtjPtr->_CSPMAyx+=yij*normx;

				PtjPtr->_CSPMAxy+=xij*normy;
				PtjPtr->_CSPMAyy+=yij*normy;
			}
		}
	}



	//2.2 ����CSPM���������Դ��
	for(unsigned int i=0;i!=Region._PtPairList.size();++i)
	{
		if(Region._PtPairList[i]._Type==enSPHPtPair)
		{
			PtiPtr=Region._PtPairList[i]._PtiPtr;
			PtjPtr=Region._PtPairList[i]._PtjPtr;
			KnlPtr=&Region._KnlList[i];

			IDi=PtiPtr->_ID;
			IDj=PtjPtr->_ID;

			mi=PtiPtr->_m;
			mj=PtjPtr->_m;

			rhoi=PtiPtr->_rho;
			rhoj=PtjPtr->_rho;

			normx=mj/rhoj*KnlPtr->_Wx;
			normy=mj/rhoj*KnlPtr->_Wy;

			bx[IDi-1]+=normx*(PtiPtr->_C-PtjPtr->_C);
			by[IDi-1]+=normy*(PtiPtr->_C-PtjPtr->_C);

			if(PtiPtr!=PtjPtr)
			{			
				normx=mi/rhoi*KnlPtr->_Wx;
				normy=mi/rhoi*KnlPtr->_Wy;

				bx[IDj-1]-=normx*(PtjPtr->_C-PtiPtr->_C);
				by[IDj-1]-=normy*(PtjPtr->_C-PtiPtr->_C);
			}
		}
	}

	//2.3 ����CSPM������ķ������򻯷���
	double mod;
	double xi;
	for(unsigned int i=0;i!=Region._PtList.size();++i)
	{
		ID=i;

		BasePtr=&Region._PtList[ID];

		ID2=BasePtr->_ID;//ID2ʵ���ϵ���i+1

		Axx=BasePtr->_CSPMAxx;
		Ayx=BasePtr->_CSPMAyx;
		Axy=BasePtr->_CSPMAxy;
		Ayy=BasePtr->_CSPMAyy;

		_OperatorBnd.Reve2ndMat(Axx,Ayx,Axy,Ayy,&Axx,&Ayx,&Axy,&Ayy);

		BasePtr->_nx=Axx*bx[ID2-1]+Ayx*by[ID2-1];
		BasePtr->_ny=Axy*bx[ID2-1]+Ayy*by[ID2-1];

		//���򻯷���
		mod=sqrt(BasePtr->_nx*BasePtr->_nx+BasePtr->_ny*BasePtr->_ny);

		xi=0.01/BasePtr->_h;

		if(mod>xi)
		{
			N[ID2-1]=1;
			BasePtr->_Nnx=BasePtr->_nx/mod;
			BasePtr->_Nny=BasePtr->_ny/mod;
		}

		else
		{
			N[ID2-1]=0;
			BasePtr->_Nnx=0.0;
			BasePtr->_Nny=0.0;
		}
	}

	bx.clear();
	by.clear();

}

void CCalBndNorm::GetWBnd(double x,double y,double h,double distance2,double & W,double & Wx,double & Wy,double & Ww)
{
	double s;
	double ttt;
	double distance;
	double dW;
	double term;
	double _norm;
	_norm=15.0/(7.0*3.1415926);

	if(distance2!=0.0)
	{
		distance=sqrt(distance2);
		s=distance/h;

		if(s<1)
		{
			W=0.66666667-s*s+0.5*s*s*s;
			dW=-2*s+1.5*s*s;
		}
		else if(s<2)
		{
			ttt=2-s;
			W=0.166666667*ttt*ttt*ttt;
			dW=-0.5*ttt*ttt;
		}
		else
		{
			W=0;
			dW=0;
		}

		W*=(_norm/(h*h));
		term=dW*_norm/(h*h*h*distance);
		Ww=term;
		Wx=term*x;
		Wy=term*y;
	}
	else
	{
		W=0.66666667*(_norm/(h*h));
		Wx=0.0;
		Wy=0.0;
	}
}
