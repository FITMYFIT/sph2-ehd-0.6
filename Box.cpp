#include"Box.h"
#include<cmath>

CBox::CBox(size_t CellNumx,size_t CellNumy)
:_CellNumx(CellNumx),_CellNumy(CellNumy),_CellList(CellNumx*CellNumy,std::vector<CBasePt *>(0,NULL))
{
}

CBox::~CBox()
{
}

void CBox::Resize(size_t CellNumx,size_t CellNumy)
{
	_CellNumx=CellNumx;
	_CellNumy=CellNumy;

	_CellList.resize(CellNumx*CellNumy,std::vector<CBasePt *>(0,NULL));
}

void CBox::UpdateBox(CRegion & Region)
{
	SetMPoint(Region);
  PutPt2Cell(Region);
}

void CBox::SetMPoint(CRegion & Region)
{
	double dp;
	unsigned int ID;

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
	dp=pow(Region._PtList[Region._CalList[0]]._Volume,0.5);//��ʱ��_CalList[0]

	_CellNumx=size_t((_MaxPointx-_MinPointx)/dp*1.0);
	_CellNumy=size_t((_MaxPointy-_MinPointy)/dp*1.0);

	if(_CellNumx>Region._ControlSPH._CellNumx)
		_CellNumx=Region._ControlSPH._CellNumx;
	if(_CellNumy>Region._ControlSPH._CellNumy)
		_CellNumy=Region._ControlSPH._CellNumy;

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
}

void CBox::PutPt2Cell(CRegion & Region)
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

void CBox::GetNbl(double x, double y, double r,std::vector<CBasePt *> & PtNbList) 
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

void CBox::GetNbl2(double x, double y, double r,std::vector<CBasePt *> & PtNbList,CRegion & Region) //�����Ա߽�����Ӷ�����
{
	double dxiac;
  double dyiac;
  double driac2;

	double dxiac2;
	double dyiac2;
	double r2;
	unsigned int ID,CellID;
	double MAXX,MINX;

	MAXX=Region._ControlSPH._PerdBndMaxX;
	MINX=Region._ControlSPH._PerdBndMinX;

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

	rsizex=size_t(floor(r/_CellSizex)+5);//���������Ѽ������ӣ���ֹ����һ���ж���պ��ӵ����
	rsizey=size_t(floor(r/_CellSizey)+1);

	for(int ii=-rsizex;ii<=rsizex;ii++)
	{
		NxI=xI+ii;
		if(NxI<0||NxI>_CellNumx)
		{
			if(NxI<0)
				NxI+=_CellNumx;
			else
				NxI-=_CellNumx;

			for(int jj=-rsizey;jj<=rsizey;jj++)
			{
				NyI=yI+jj;
				if(NyI>=0&&NyI<_CellNumy)
				{
					CellID=NyI*_CellNumx+NxI;

					for(int ff=0;ff<_CellList[NyI*_CellNumx+NxI].size();ff++)
					{
						if(_CellList[CellID][ff]->_x>x)//�������������߽�
							dxiac=_CellList[CellID][ff]->_x-(MAXX-MINX)-x;						
						if(_CellList[CellID][ff]->_x<x)//�����������Ҳ�߽�
							dxiac=_CellList[CellID][ff]->_x+(MAXX-MINX)-x;

						dyiac=_CellList[CellID][ff]->_y-y;

						dxiac2=dxiac*dxiac;
						dyiac2=dyiac*dyiac;

						driac2=dxiac2+dyiac2;

						if(driac2<=r2)
						{ 
							PtNbList.push_back(_CellList[CellID][ff]);
						}
					}	
				}
			}
		}
	}}

void CBox::GetNearBndPt(std::vector<CBasePt *> & BndPt)//�ҳ��߽總��������
{
	unsigned int i,j,k;
	vector<double> RowMinX;//ÿһ������x�������Сֵ�����������Ա߽�����
	vector<double> RowMaxX;//ÿһ������x��������ֵ�����������Ա߽�����
	vector<int> RowEmptyId;//����ָʾÿһ���ǲ���ȫΪ�պ���
	unsigned int NBCELLX,NBCELLY;

	NBCELLX=_CellNumx;
	NBCELLY=_CellNumy;


	RowMinX.resize(NBCELLY,0.0);
	RowMaxX.resize(NBCELLY,0.0);
	RowEmptyId.resize(NBCELLY,0);

	//�����Ա߽�1.�ҳ�ÿһ�е������������ֵ����Сֵ
	for(i=0;i<NBCELLY;i++)//iΪ�к�
	{
		//Ϊÿһ�е������������С�����ֵ����ֵ
		for(j=0;j<NBCELLX;j++)
		{
			if(_CellList[i*NBCELLX+j].size()!=0)
			{
				RowMinX[i]=_CellList[i*NBCELLX+j][0]->_x;
				RowMaxX[i]=_CellList[i*NBCELLX+j][0]->_x;

				RowEmptyId[i]=1;//������һ�в�ȫ�ǿպ���

				break;
			}
		}

		if(RowEmptyId[i]==1)
		{
			for(j=0;j<NBCELLX;j++)
			{
				for(k=1;k<_CellList[i*NBCELLX+j].size();k++)
				{
					if(_CellList[i*NBCELLX+j][k]->_x<RowMinX[i])
						RowMinX[i]=_CellList[i*NBCELLX+j][k]->_x;
					if(_CellList[i*NBCELLX+j][k]->_x>RowMinX[i])
						RowMaxX[i]=_CellList[i*NBCELLX+j][k]->_x;
				}
			}
		}
	}

	//�����Ա߽�2.�ж���Щ�����Ǵ��ڱ߽總��������
	
	BndPt.clear();//�������

	unsigned int nbbndpt=0;
	unsigned int kk=0;
	for(i=0;i<NBCELLY;i++)
	{
		for(j=0;j<NBCELLX;j++)
		{
			for(k=0;k<_CellList[i*NBCELLX+j].size();k++)
			{
				if(_CellList[i*NBCELLX+j][k]->_x-RowMinX[i]<_CellList[i*NBCELLX+j][k]->_r)
				{
					BndPt.push_back(_CellList[i*NBCELLX+j][k]);
				}		

				if(RowMaxX[i]-_CellList[i*NBCELLX+j][k]->_x<_CellList[i*NBCELLX+j][k]->_r)
				{
					BndPt.push_back(_CellList[i*NBCELLX+j][k]);
				}
			}
		}
	}
}
