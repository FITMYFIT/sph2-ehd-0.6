#include "CalculateRange.h"

CCalculateRange::CCalculateRange()
{
}

CCalculateRange::~CCalculateRange()
{
}


//˵����������������������ӵ�_Iflag=1;ֻ�����������������������±�ţ��Խ�����ʽ���㣩��Region._CalListֻ������������������ӵı��
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

	//���������ӵ�ID2��Iflag��0
	for(unsigned int i=0;i!=Region._PtList.size();++i)
	{
		Region._PtList[i]._ID2=0;
		Region._PtList[i]._Iflag=0;
	}

	Region._CalList.clear();

	//1. X Y�������򶼲�����
	if((Xcrmin==0.0&&Xcrmax==0.0)&&(Ycrmin==0.0&&Ycrmax==0.0))
	{
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			if(Region._PtList[i]._Type==enSPHPt)//�������ӣ�������ָʾ��1�����±�ţ��Խ�����ʽ���㣩
			{
				icount++;
				Region._PtList[i]._ID2=icount;
				Region._PtList[i]._Iflag=1;

				Region._CalList.push_back(i);
			}

			//if(Region._PtList[i]._Type==enNULLPt||Region._PtList[i]._Type==enBndPt||Region._PtList[i]._Type==enDumPt)//Null���ӣ�������ָʾ��1�������±�ţ���������ʽ���㣩
      else
        {
          Region._PtList[i]._ID2=0;
          Region._PtList[i]._Iflag=1;
        }
		}
	}

	//2. X�������ޣ�Y��������
	else if((Xcrmin==0.0&&Xcrmax==0.0)&&((Ycrmin==0.0&&Ycrmax==0.0)==0))
	{
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			BasePtPtr=&Region._PtList[i];

			if(BasePtPtr->_y>=Ycrmin&&BasePtPtr->_y<=Ycrmax)
			{
				if(BasePtPtr->_Type==enSPHPt)//�������ӣ����±�ţ��Խ�����ʽ���㣩��������ָʾ��1
				{
					icount++;
					BasePtPtr->_ID2=icount;//�����Ӱ��Ƿ������������±�ţ����ڽ�����ʽ����
					BasePtPtr->_Iflag=1;

					Region._CalList.push_back(i);
				}

        // if(Region._PtList[i]._Type==enNULLPt||Region._PtList[i]._Type==enBndPt||Region._PtList[i]._Type==enDumPt)//Null���ӣ�������ָʾ��1�������±�ţ���������ʽ���㣩
        else
            {
              BasePtPtr->_ID2=0;//�����Ӱ��Ƿ������������±�ţ����ڽ�����ʽ����
              BasePtPtr->_Iflag=1;
            }
			}
		}
	}

	//3. Y�������ޣ�X��������
	else if((Ycrmin==0.0&&Ycrmax==0.0)&&((Xcrmin==0.0&&Xcrmax==0.0)==0))
	{
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			BasePtPtr=&Region._PtList[i];

			if(BasePtPtr->_x>=Xcrmin&&BasePtPtr->_x<=Xcrmax)
			{
				if(BasePtPtr->_Type==enSPHPt)//�������ӣ����±�ţ��Խ�����ʽ���㣩��������ָʾ��1
				{
					icount++;
					BasePtPtr->_ID2=icount;//�����Ӱ��Ƿ������������±�ţ����ڽ�����ʽ����
					BasePtPtr->_Iflag=1;

					Region._CalList.push_back(i);
				}

        // if(Region._PtList[i]._Type==enNULLPt||Region._PtList[i]._Type==enBndPt||Region._PtList[i]._Type==enDumPt)//Null���ӣ�������ָʾ��1�������±�ţ���������ʽ���㣩
        else
          {
            BasePtPtr->_ID2=0;//�����Ӱ��Ƿ������������±�ţ����ڽ�����ʽ����
            BasePtPtr->_Iflag=1;
          }
			}
		}
	}

	//4. X Y������������
	else
	{
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			BasePtPtr=&Region._PtList[i];

			if(BasePtPtr->_x>=Xcrmin&&BasePtPtr->_x<=Xcrmax
				&&BasePtPtr->_y>=Ycrmin&&BasePtPtr->_y<=Ycrmax)
			{
				if(BasePtPtr->_Type==enSPHPt)//�������ӣ����±�ţ��Խ�����ʽ���㣩��������ָʾ��1
				{
					icount++;
					BasePtPtr->_ID2=icount;//�����Ӱ��Ƿ������������±�ţ����ڽ�����ʽ����
					BasePtPtr->_Iflag=1;

					Region._CalList.push_back(i);
				}

        //  if(Region._PtList[i]._Type==enNULLPt||Region._PtList[i]._Type==enBndPt||Region._PtList[i]._Type==enDumPt)//Null���ӣ�������ָʾ��1�������±�ţ���������ʽ���㣩
        else
          {
            BasePtPtr->_ID2=0;//�����Ӱ��Ƿ������������±�ţ����ڽ�����ʽ����
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
