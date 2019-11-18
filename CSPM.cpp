#include "CSPM.h"

CCSPM::CCSPM()
{
}

CCSPM::~CCSPM()
{
}

void CCSPM::GetCSPMFunCorctCoef(CRegion &Region)
{
	unsigned int ID;
	CBasePt * PtiPtr,*PtjPtr;
	CKnl * KnlPtr;

	//2.求出CSPM修正系数 sum(mj/rhoj*wij)
	for(unsigned int i=0;i!=Region._PtPairList.size();++i)
	{
		PtiPtr=Region._PtPairList[i]._PtiPtr;
		PtjPtr=Region._PtPairList[i]._PtjPtr;
		KnlPtr=&Region._KnlList[i];

		if(Region._PtPairList[i]._Type==enSPHPtPair||Region._PtPairList[i]._Type==enSPHNULLPtPair)
		{
			PtiPtr->_CSPMCoef+=PtjPtr->_m/PtjPtr->_rho*KnlPtr->_W;

			if(PtiPtr!=PtjPtr)
			{
				PtjPtr->_CSPMCoef+=PtiPtr->_m/PtiPtr->_rho*KnlPtr->_W;
			}
		}
	}

	Region._ControlSPH._CSPMIflag1=1;//监控参数置1，表示已经计算过了
}

//Ref: SPH_4 P119
//one more step: calculate the inverse of matrix A, so can be used directly multiplies b, 2019.10.11
void CCSPM::GetCSPMGradCorctCoef(CRegion &Region)//计算CSPM对函数导数修正时的系数矩阵
{
	//这里系数矩阵的表示方法与正常不太一样，（主要指x y z 的写法），这是与CSPM具体表达式中的写法一致的
	//[ Axx Ayx Azx ]   [fx]   [bx]
	//[ Axy Ayy Azy ] * [fy] = [by]
	//[ Axz Ayz Azz ]   [fz]   [bz]
	// Axx=sum(mj/rhoj*xij*gradwx)
	// Ayx=sum(mj/rhoj*yij*gradwx)
	// Azx=sum(mj/rhoj*zij*gradwx)
	// Axy=sum(mj/rhoj*xij*gradwy)
	// Ayy=sum(mj/rhoj*yij*gradwy)
	// Azy=sum(mj/rhoj*zij*gradwy)
	// Axz=sum(mj/rhoj*xij*gradwz)
	// Ayz=sum(mj/rhoj*yij*gradwz)
	// Azz=sum(mj/rhoj*zij*gradwz)


	CBasePt * PtiPtr,* PtjPtr;
  CBasePt * PtPtr;
	CKnl * KnlPtr;

	double norm;
	double mi,mj;
	double rhoi,rhoj;
	double xij,yij;
	double normx,normy;

  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      PtPtr=&Region._PtList[i];
      PtPtr->_CSPMAxx=0.0;
      PtPtr->_CSPMAyx=0.0;
      PtPtr->_CSPMAzx=0.0;
      PtPtr->_CSPMAxy=0.0;
      PtPtr->_CSPMAyy=0.0;
      PtPtr->_CSPMAzy=0.0;
      PtPtr->_CSPMAxz=0.0;
      PtPtr->_CSPMAyz=0.0;
      PtPtr->_CSPMAzz=0.0;
    }

  for (size_t i=0;i<Region._PtPairList.size();i++)
      {
        if(Region._PtPairList[i]._Type==enSPHPtPair||Region._PtPairList[i]._Type==enSPHNULLPtPair)
          {
            PtiPtr=Region._PtPairList[i]._PtiPtr;
            PtjPtr=Region._PtPairList[i]._PtjPtr;
            KnlPtr=&Region._KnlList[i];

            mi=PtiPtr->_m;
            mj=PtjPtr->_m;

            rhoi=PtiPtr->_rho;
            rhoj=PtjPtr->_rho;

            xij=Region._PtPairList[i]._xij;
            yij=Region._PtPairList[i]._yij;

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

    //calculate the inverse of coefficient matrix A 2019.10.11
    CBasePt * BasePtPtr;
    COperator Operator;
    double Axx,Axy,Ayx,Ayy;
    for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
 
      Axx=BasePtPtr->_CSPMAxx;
      Axy=BasePtPtr->_CSPMAxy;
      Ayx=BasePtPtr->_CSPMAyx;
      Ayy=BasePtPtr->_CSPMAyy;
      Operator.Reve2ndMat(Axx,Ayx,Axy,Ayy,
                          &BasePtPtr->_CSPMAxx,&BasePtPtr->_CSPMAyx,&BasePtPtr->_CSPMAxy,&BasePtPtr->_CSPMAyy);
    }

  Region._ControlSPH._CSPMIflag2=1;//监控参数置1，表示已经计算过了
}

//Ref: SPH_4 P119
//one more step: calculate the inverse of matrix A, so can be used directly multiplies b, 2019.10.11
void CCSPM::GetCSPMGradCorctCoef2(CRegion &Region)//计算CSPM对函数导数修正时的系数矩阵
{
	//这里系数矩阵的表示方法与正常不太一样，（主要指x y z 的写法），这是与CSPM具体表达式中的写法一致的
	//[ Axx Ayx Azx ]   [fx]   [bx]
	//[ Axy Ayy Azy ] * [fy] = [by]
	//[ Axz Ayz Azz ]   [fz]   [bz]
	// Axx=sum(mj/rhoj*xij*gradwx)
	// Ayx=sum(mj/rhoj*yij*gradwx)
	// Azx=sum(mj/rhoj*zij*gradwx)
	// Axy=sum(mj/rhoj*xij*gradwy)
	// Ayy=sum(mj/rhoj*yij*gradwy)
	// Azy=sum(mj/rhoj*zij*gradwy)
	// Axz=sum(mj/rhoj*xij*gradwz)
	// Ayz=sum(mj/rhoj*yij*gradwz)
	// Azz=sum(mj/rhoj*zij*gradwz)


	CBasePt * PtiPtr,* PtjPtr;
  CBasePt * PtPtr;
	CKnl * KnlPtr;

	double norm;
	double mi,mj;
	double rhoi,rhoj;
	double xij,yij;
	double normx,normy;

  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      PtPtr=&Region._PtList[i];
      PtPtr->_CSPMAxx=0.0;
      PtPtr->_CSPMAyx=0.0;
      PtPtr->_CSPMAzx=0.0;
      PtPtr->_CSPMAxy=0.0;
      PtPtr->_CSPMAyy=0.0;
      PtPtr->_CSPMAzy=0.0;
      PtPtr->_CSPMAxz=0.0;
      PtPtr->_CSPMAyz=0.0;
      PtPtr->_CSPMAzz=0.0;
    }

  for (size_t i=0;i<Region._PtPairList.size();i++)
      {
        if(Region._PtPairList[i]._Type==enSPHPtPair
           ||Region._PtPairList[i]._Type==enSPHNULLPtPair
           ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
           ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
          {
            PtiPtr=Region._PtPairList[i]._PtiPtr;
            PtjPtr=Region._PtPairList[i]._PtjPtr;
            KnlPtr=&Region._KnlList[i];

            mi=PtiPtr->_m;
            mj=PtjPtr->_m;

            rhoi=PtiPtr->_rho;
            rhoj=PtjPtr->_rho;

            xij=Region._PtPairList[i]._xij;
            yij=Region._PtPairList[i]._yij;

            normx=mj/rhoj*KnlPtr->_Wx;
            normy=mj/rhoj*KnlPtr->_Wy;

            PtiPtr->_CSPMAxx+=xij*normx;
            PtiPtr->_CSPMAyx+=yij*normx;
		
            PtiPtr->_CSPMAxy+=xij*normy;
            PtiPtr->_CSPMAyy+=yij*normy;

            if(PtiPtr!=PtjPtr&&PtjPtr->_Type==enSPHPt)
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

  //calculate the inverse of coefficient matrix A 2019.10.11
  CBasePt * BasePtPtr;
  COperator Operator;
  double Axx,Axy,Ayx,Ayy;
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
 
      Axx=BasePtPtr->_CSPMAxx;
      Axy=BasePtPtr->_CSPMAxy;
      Ayx=BasePtPtr->_CSPMAyx;
      Ayy=BasePtPtr->_CSPMAyy;
      Operator.Reve2ndMat(Axx,Ayx,Axy,Ayy,
                          &BasePtPtr->_CSPMAxx,&BasePtPtr->_CSPMAyx,&BasePtPtr->_CSPMAxy,&BasePtPtr->_CSPMAyy);
    }

  Region._ControlSPH._CSPMIflag2=2;//监控参数置1，表示已经计算过了
}
