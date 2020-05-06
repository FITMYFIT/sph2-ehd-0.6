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

	//找到第一个处于计算域中的粒子（流体粒子或Null粒子）
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

	//找到计算域内的最小、最大值
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

	_MinPointr=Region._PtList[0]._r;//暂时的，肯定有问题

	_MaxPointx+=_MinPointr*0.5;
	_MinPointx-=_MinPointr*0.5;

	_MaxPointy+=_MinPointr*0.5;
	_MinPointy-=_MinPointr*0.5;

	//动态控制每步布的盒子数
	dp=pow(Region._PtList[Region._CalList[0]]._Volume,0.5);//暂时用_CalList[0]

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
	
	//重置盒子数并将内部清空
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

void CBox::GetNbl2(double x, double y, double r,std::vector<CBasePt *> & PtNbList,CRegion & Region) //周期性边界的粒子对搜索
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

	rsizex=size_t(floor(r/_CellSizex)+5);//比正常多搜几个盒子，防止出现一行有多个空盒子的情况
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
						if(_CellList[CellID][ff]->_x>x)//主搜粒子在左侧边界
							dxiac=_CellList[CellID][ff]->_x-(MAXX-MINX)-x;						
						if(_CellList[CellID][ff]->_x<x)//主搜粒子在右侧边界
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

void CBox::GetNearBndPt(std::vector<CBasePt *> & BndPt)//找出边界附近的粒子
{
	unsigned int i,j,k;
	vector<double> RowMinX;//每一行粒子x坐标的最小值，用于周期性边界条件
	vector<double> RowMaxX;//每一行粒子x坐标的最大值，用于周期性边界条件
	vector<int> RowEmptyId;//用以指示每一行是不是全为空盒子
	unsigned int NBCELLX,NBCELLY;

	NBCELLX=_CellNumx;
	NBCELLY=_CellNumy;


	RowMinX.resize(NBCELLY,0.0);
	RowMaxX.resize(NBCELLY,0.0);
	RowEmptyId.resize(NBCELLY,0);

	//周期性边界1.找出每一行的粒子坐标最大值和最小值
	for(i=0;i<NBCELLY;i++)//i为行号
	{
		//为每一行的粒子坐标的最小、最大值赋初值
		for(j=0;j<NBCELLX;j++)
		{
			if(_CellList[i*NBCELLX+j].size()!=0)
			{
				RowMinX[i]=_CellList[i*NBCELLX+j][0]->_x;
				RowMaxX[i]=_CellList[i*NBCELLX+j][0]->_x;

				RowEmptyId[i]=1;//表明这一行不全是空盒子

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

	//周期性边界2.判断哪些粒子是处于边界附近的粒子
	
	BndPt.clear();//清空容器

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
