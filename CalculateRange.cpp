#include "CalculateRange.h"

CCalculateRange::CCalculateRange()
{
}

CCalculateRange::~CCalculateRange()
{
}


//说明：令进入计算域的所有粒子的_Iflag=1;只将进入计算域的流体粒子重新编号（以进行隐式计算），Region._CalList只存进入计算域的流体粒子的编号
void CCalculateRange::_UpdateRange(CRegion& Region,unsigned int TimeSteps)
{
	double Xcrmin,Xcrmax,Ycrmin,Ycrmax;

	CBasePt * BasePtPtr;
	CPart * PartPtr;
	CSPHPt * PtPtr;

	unsigned int icount=0;

	Xcrmin=Region._ControlSPH._Xcrmin;
	Xcrmax=Region._ControlSPH._Xcrmax;
	Ycrmin=Region._ControlSPH._Ycrmin;
	Ycrmax=Region._ControlSPH._Ycrmax;

	//将所有粒子的ID2和Iflag置0
	for(unsigned int i=0;i!=Region._PtList.size();++i)
	{
		Region._PtList[i]._ID2=0;
		Region._PtList[i]._Iflag=0;
	}

	Region._CalList.clear();

	//1. X Y两个方向都不设限
	if((Xcrmin==0.0&&Xcrmax==0.0)&&(Ycrmin==0.0&&Ycrmax==0.0))
	{
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			if(Region._PtList[i]._Type==enSPHPt)//流体粒子，计算域指示置1，重新编号（以进行隐式计算）
			{
				icount++;
				Region._PtList[i]._ID2=icount;
				Region._PtList[i]._Iflag=1;

				Region._CalList.push_back(i);
			}

			//if(Region._PtList[i]._Type==enNULLPt||Region._PtList[i]._Type==enBndPt||Region._PtList[i]._Type==enDumPt)//Null粒子，计算域指示置1，不重新编号（不进行隐式计算）
      else
        {
          Region._PtList[i]._ID2=0;
          Region._PtList[i]._Iflag=1;
        }
		}
	}

	//2. X方向不设限，Y方向设限
	else if((Xcrmin==0.0&&Xcrmax==0.0)&&((Ycrmin==0.0&&Ycrmax==0.0)==0))
	{
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			BasePtPtr=&Region._PtList[i];

			if(BasePtPtr->_y>=Ycrmin&&BasePtPtr->_y<=Ycrmax)
			{
				if(BasePtPtr->_Type==enSPHPt)//流体粒子，重新编号（以进行隐式计算），计算域指示置1
				{
					icount++;
					BasePtPtr->_ID2=icount;//将粒子按是否进入计算域重新编号，便于进行隐式计算
					BasePtPtr->_Iflag=1;

					Region._CalList.push_back(i);
				}

        // if(Region._PtList[i]._Type==enNULLPt||Region._PtList[i]._Type==enBndPt||Region._PtList[i]._Type==enDumPt)//Null粒子，计算域指示置1，不重新编号（不进行隐式计算）
        else
            {
              BasePtPtr->_ID2=0;//将粒子按是否进入计算域重新编号，便于进行隐式计算
              BasePtPtr->_Iflag=1;
            }
			}
		}
	}

	//3. Y方向不设限，X方向设限
	else if((Ycrmin==0.0&&Ycrmax==0.0)&&((Xcrmin==0.0&&Xcrmax==0.0)==0))
	{
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			BasePtPtr=&Region._PtList[i];

			if(BasePtPtr->_x>=Xcrmin&&BasePtPtr->_x<=Xcrmax)
			{
				if(BasePtPtr->_Type==enSPHPt)//流体粒子，重新编号（以进行隐式计算），计算域指示置1
				{
					icount++;
					BasePtPtr->_ID2=icount;//将粒子按是否进入计算域重新编号，便于进行隐式计算
					BasePtPtr->_Iflag=1;

					Region._CalList.push_back(i);
				}

        // if(Region._PtList[i]._Type==enNULLPt||Region._PtList[i]._Type==enBndPt||Region._PtList[i]._Type==enDumPt)//Null粒子，计算域指示置1，不重新编号（不进行隐式计算）
        else
          {
            BasePtPtr->_ID2=0;//将粒子按是否进入计算域重新编号，便于进行隐式计算
            BasePtPtr->_Iflag=1;
          }
			}
		}
	}

	//4. X Y两个方向都设限
	else
	{
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			BasePtPtr=&Region._PtList[i];

			if(BasePtPtr->_x>=Xcrmin&&BasePtPtr->_x<=Xcrmax
				&&BasePtPtr->_y>=Ycrmin&&BasePtPtr->_y<=Ycrmax)
			{
				if(BasePtPtr->_Type==enSPHPt)//流体粒子，重新编号（以进行隐式计算），计算域指示置1
				{
					icount++;
					BasePtPtr->_ID2=icount;//将粒子按是否进入计算域重新编号，便于进行隐式计算
					BasePtPtr->_Iflag=1;

					Region._CalList.push_back(i);
				}

        //  if(Region._PtList[i]._Type==enNULLPt||Region._PtList[i]._Type==enBndPt||Region._PtList[i]._Type==enDumPt)//Null粒子，计算域指示置1，不重新编号（不进行隐式计算）
        else
          {
            BasePtPtr->_ID2=0;//将粒子按是否进入计算域重新编号，便于进行隐式计算
            BasePtPtr->_Iflag=1;
          }
			}
		}
	}

	Region._StatDataList._InvolvedFluidNum=Region._CalList.size();
}

void CCalculateRange::_Clear(CRegion &Region)
{
	Region._CalList.clear();

	//for(unsigned int i=0;i!=Region._PartList.size();++i)
	//{
	//	Region._PartList[i]._PartCalList.clear();
	//}
}
