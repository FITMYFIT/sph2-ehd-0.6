#include"GetKnlList.h"

CGetKnlList::CGetKnlList(enKNLTYPE Type)
:_Type(Type)
{
	_norm=15.0/(7.0*3.1415926);
}

CGetKnlList::~CGetKnlList()
{
}

void CGetKnlList::GetKnlList(CRegion &Region)
{
	CKnl *KnlPtr;
	CKnl *Knl2Ptr;
	CBasePt * PtiPtr,* PtjPtr;
	double x,y,h,H,distance2;
	unsigned int icount=0;
	double MAXX,MINX;//用于周期性边界条件
	double DisCrit;//用于周期性边界条件，判断两粒子是否是利用周期性边界相互作用

	if (Region._ControlSPH._PerdBnd==1)
	{	
		MAXX=Region._ControlSPH._PerdBndMaxX;
		MINX=Region._ControlSPH._PerdBndMinX;
		DisCrit=Region._PtList[0]._r*2;
	}

	Region._KnlList.resize(Region._PtPairList.size());

	icount=0;
	for(size_t i=0;i<Region._PtPairList.size();i++)
	{
		//if(Region._PtPairList[i]._Type==enSPHPtPair||Region._PtPairList[i]._Type==enSPHNULLPtPair
    //      ||Region._PtPairList[i]._Type==enSPHBndPtPair||Region._PtPairList[i]._Type==enSPHDumPtPair)
		{
			PtiPtr=Region._PtPairList[i]._PtiPtr;
			PtjPtr=Region._PtPairList[i]._PtjPtr;

      x=Region._PtPairList[i]._xij;
      y=Region._PtPairList[i]._yij;	
      h=0.5*(PtiPtr->_h+PtjPtr->_h); 

      if (Region._ControlSPH._PerdBnd==1)//周期性边界条件
          {		
            if(x>DisCrit)
              x=x-(MAXX-MINX);
            if(x<-DisCrit)
              x=x+(MAXX-MINX);
          }

        distance2=Region._PtPairList[i]._driac2;
			
			GetW(x,y,h,distance2,Region._KnlList[icount]._W,Region._KnlList[icount]._Wx,Region._KnlList[icount]._Wy,Region._KnlList[icount]._Ww);

			////1.用粒子j的光滑长度
			//h=PtjPtr->_h;
			//GetW(x,y,h,distance2,Region._KnlList[icount]._W,Region._KnlList[icount]._Wx,Region._KnlList[icount]._Wy,Region._KnlList[icount]._Ww);
			////2.用粒子i的光滑长度
			//h=PtiPtr->_h;
			//GetW(x,y,h,distance2,Region._KnlList[icount]._Wj,Region._KnlList[icount]._Wxj,Region._KnlList[icount]._Wyj,Region._KnlList[icount]._Wwj);

			icount++;
		}		
	}

	if(Region._ControlSPH._SPHST==1)
	{
		icount=0;
		Region._KnlList2.resize(Region._PtPairList2.size());

		for(size_t i=0;i<Region._PtPairList2.size();i++)
		{
			if(Region._PtPairList2[i]._Type==enSPHPtPair||Region._PtPairList2[i]._Type==enSPHNULLPtPair||Region._PtPairList2[i]._Type==enSPHBndPtPair)
			{
				PtiPtr=(CSPHPt*)Region._PtPairList2[i]._PtiPtr;
				PtjPtr=(CSPHPt*)Region._PtPairList2[i]._PtjPtr;

				x=PtiPtr->_x-PtjPtr->_x;
				y=PtiPtr->_y-PtjPtr->_y;
				H=0.5*(PtiPtr->_H+PtjPtr->_H);
				distance2=Region._PtPairList2[i]._driac2;

				GetW(x,y,H,distance2,Region._KnlList2[icount]._W,Region._KnlList2[icount]._Wx,Region._KnlList2[icount]._Wy,Region._KnlList2[icount]._Ww);
				icount++;
			}		
		}	
	}
}

void CGetKnlList::GetW(double x,double y,double h,double distance2,double & W,double & Wx,double & Wy,double & Ww)
{
	double s;
	double ttt;
	double distance;
	double dW;
	double term;

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

void CGetKnlList::GetW2(double x,double y,double h,double distance2,double & W,double & Wx,double & Wy,double & W2)
{
	double s;
	double ttt;
	double distance;
	double dW;
	double term;

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

void CGetKnlList::GetW0(double distance0,double h,double & W)
{
	double s;
	double ttt;

	if(distance0!=0.0)
	{
		s=distance0/h;

		if(s<1)
		{
			W=0.66666667-s*s+0.5*s*s*s;
		}
		else if(s<2)
		{
			ttt=2-s;
			W=0.166666667*ttt*ttt*ttt;
		}
		else
		{
			W=0;
		}
		W*=(_norm/(h*h));
	}
	else
	{
		W=0.66666667*(_norm/(h*h));
	}
}


void CGetKnlList::ClearKnlList(CRegion &Region)
{
	////delete KnlPtr;
	//for(size_t j=0;j<Region._KnlList.size();j++)
	//{
	//	if(Region._KnlList[j]!=NULL)
	//	{
	//		delete Region._KnlList[j];
	//		Region._KnlList[j]=NULL;
	//	}
	//}

	//for(size_t j=0;j<Region._KnlList2.size();j++)
	//{
	//	if(Region._KnlList2[j]!=NULL)
	//	{
	//		delete Region._KnlList2[j];
	//		Region._KnlList2[j]=NULL;
	//	}
	//}

	Region._KnlList.clear();

	if(Region._ControlSPH._SPHST==1)
	{
		Region._KnlList2.clear();
	}
}


void CGetKnlList::GetCrctKnlGrad(CRegion & Region,CBasePt * PtPtr,CKnl * KnlPtr,double & wcx,double & wcy)//修正核函数梯度
{
	if(Region._ControlSPH._CSPMIflag2==0)//如果没有计算CSPM修正参数，则需要计算
	{
		_CSPMG.GetCSPMGradCorctCoef(Region);
	}

	double L11,L12,L21,L22;

	L11=-PtPtr->_CSPMAxx;   L12=-PtPtr->_CSPMAyx;
	L21=-PtPtr->_CSPMAxy;   L22=-PtPtr->_CSPMAyy;

	_OperatorG.Reve2ndMat(L11,L12,L21,L22,
											&L11,&L12,&L21,&L22);

	wcx=L11*KnlPtr->_Wx+L12*KnlPtr->_Wy;
	wcy=L21*KnlPtr->_Wx+L22*KnlPtr->_Wy;
}

void CGetKnlList::GetMshKnlList( CRegion &Region )
{
	CKnl *KnlPtr;
	CKnl *Knl2Ptr;
	CBasePt * PtiPtr,* PtjPtr;
	CMesh * MshjPtr;
	double x,y,h,H,distance2;
	unsigned int icount=0;
	double MAXX,MINX;//用于周期性边界条件

	if (Region._ControlSPH._PerdBnd==1)
	{	
		MAXX=Region._ControlSPH._PerdBndMaxX;
		MINX=Region._ControlSPH._PerdBndMinX;
	}

	Region._PtMshKnlList.resize(Region._PtMshPairList.size());

	icount=0;
	for(size_t i=0;i<Region._PtMshPairList.size();i++)
	{
		PtiPtr=Region._PtMshPairList[i]._PtiPtr;
		MshjPtr=Region._PtMshPairList[i]._MshjPtr;

		x=PtiPtr->_x-MshjPtr->_x;
		y=PtiPtr->_y-MshjPtr->_y;
		//h=0.5*(PtiPtr->_h+MshjPtr->_h); 
		h=PtiPtr->_h;//用粒子的光滑长度

		if (Region._ControlSPH._PerdBnd==1)//周期性边界条件
		{		
			if(x>4*h)
				x=x-(MAXX-MINX);
			if(x<-4*h)
				x=x+(MAXX-MINX);
		}

		distance2=Region._PtMshPairList[i]._driac2;

		GetW(x,y,h,distance2,Region._PtMshKnlList[icount]._W,Region._PtMshKnlList[icount]._Wx,Region._PtMshKnlList[icount]._Wy,Region._PtMshKnlList[icount]._Ww);
		icount++;
	}
}
