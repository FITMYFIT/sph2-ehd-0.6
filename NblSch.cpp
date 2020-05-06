#include "NblSch.h"
#include<cmath>

CNblSch::CNblSch(size_t CellNumx,size_t CellNumy)
:_Box(CellNumx,CellNumy),_PtNbList(0,NULL)
{
}

CNblSch::~CNblSch()
{
}

void CNblSch::ResizeBox(size_t CellNumx,size_t CellNumy/*,size_t CellNumz*/)
{
	_Box.Resize(CellNumx,CellNumy/*,CellNumz*/);
}


 void CNblSch::GetNbl(CRegion& Region)
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
  CBasePt * PtPtr;
  double MAXX,MINX;
	//if (Region._ControlSPH._IFlagRemesh==0)
	{
		_Box.UpdateBox(Region);//update the Box
	}

  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      Region._PtList[i]._NumNegbor=0;
    }

  //周期性边界，目前只能用于X方向的周期性边界
  //周期性边界1.找出每一行的粒子坐标最大值和最小值
  if(Region._ControlSPH._PerdBnd==1)
    {
      MAXX=Region._ControlSPH._PerdBndMaxX;
      MINX=Region._ControlSPH._PerdBndMinX;

      ////找出边界附近的粒子
      //_Box.GetNearBndPt(BndPt);

      //周期性边界3.周期性边界的粒子对搜索
      for(unsigned int f=0;f<Region._PtList.size();f++)
        {
          PtPtr=&Region._PtList[f];
          if(PtPtr->_Type==enSPHPt)
            {
              _Box.GetNbl2(PtPtr->_x,PtPtr->_y,PtPtr->_r,_PtNbList,Region);
				
              for(unsigned int i=0;i!=_PtNbList.size();++i)
                {
                  if(_PtNbList[i]->_x>PtPtr->_x)//主粒子在左侧边界附近
                    dxiac=PtPtr->_x-_PtNbList[i]->_x+(MAXX-MINX);
                  else//主粒子在右侧边界附近
                    dxiac=PtPtr->_x-_PtNbList[i]->_x-(MAXX-MINX);

                  dyiac=PtPtr->_y-_PtNbList[i]->_y;
                  dxiac2=dxiac*dxiac;
                  dyiac2=dyiac*dyiac;
                  driac2=dxiac2+dyiac2;
                  r=_PtNbList[i]->_r+PtPtr->_r;
                  r2=r*r;

                  //the second searching
                  if(4*driac2<=r2)
                    {
                      switch (_PtNbList[i]->_Type)
                        {
                        case enSPHPt: 	
                          if(Region._PtList[f]._ID<=_PtNbList[i]->_ID)
                            {
                              PushPtPair(Region._PtPairList, PtPtr, _PtNbList[i], enSPHPtPair, dxiac, dyiac, driac2);
                            }
                          //if the PartID of the pt which is searching is smaller or the same as  the searched one ,put it inside the _PtPairList	
                          else
                            {
                              double PtNbr;
                              PtNbr=_PtNbList[i]->_r*_PtNbList[i]->_r;
                              if(driac2>PtNbr)
                                {
                                  PushPtPair(Region._PtPairList, _PtNbList[i], PtPtr, enSPHPtPair, -dxiac, -dyiac, driac2);
                                }
                            } 
                          break;
                          
                        case enNULLPt:
                          PushPtPair(Region._PtPairList, PtPtr, _PtNbList[i], enSPHNULLPtPair, dxiac, dyiac, driac2);
                          break;

                        case enBndPt:
                          PushPtPair(Region._PtPairList, PtPtr, _PtNbList[i], enSPHBndPtPair, dxiac, dyiac, driac2);
                          break;

                        case enGhost1Pt:
                          PushPtPair(Region._PtPairList, PtPtr, _PtNbList[i], enSPHGhostPtPair, dxiac, dyiac, driac2);
                          break;

                        case enEHDBndPt://2019.10.26
                          PushPtPair(Region._PtPairList, PtPtr, _PtNbList[i], enSPHEHDBndPtPair, dxiac, dyiac, driac2);
                          break;
                    
                        case enEHDDumPt://2019.10.26
                          PushPtPair(Region._PtPairList, PtPtr, _PtNbList[i], enSPHEHDDumPtPair, dxiac, dyiac, driac2);
                          break;

                        default:
                          break;
                        }
                    }
                }
            }
        }
    }

	//正常找粒子对部分
	//只用流体粒子作为主搜粒子
	for(size_t f=0;f<Region._CalList.size();f++)
    {
      ID=Region._CalList[f];

      PtPtr=&Region._PtList[ID];
      
      if(PtPtr->_Type==enSPHPt)
        {
          vector<CBasePt *> PtNbList;
          _Box.GetNbl(PtPtr->_x,PtPtr->_y,PtPtr->_r, PtNbList);

          for(size_t i=0;i<PtNbList.size();i++)//for all the number of the PtNbList
            {
              dxiac=PtPtr->_x-PtNbList[i]->_x;
              dyiac=PtPtr->_y-PtNbList[i]->_y;

              dxiac2=dxiac*dxiac;
              dyiac2=dyiac*dyiac;

              driac2=dxiac2+dyiac2;
              r=PtNbList[i]->_r+PtPtr->_r;
              r2=r*r;

              //the second searching
              if(4*driac2<=r2)
                {						
                  switch (PtNbList[i]->_Type)
                    {
                    case enSPHPt: 
                      {	
                        if((ID+1)<=PtNbList[i]->_ID)
                          {
                            PushPtPair(Region._PtPairList, PtPtr, PtNbList[i], enSPHPtPair, dxiac, dyiac, driac2);   
                          }

                        else
                          {
                            double PtNbr;
                            PtNbr=PtNbList[i]->_r*PtNbList[i]->_r;
                            if(driac2>PtNbr)
                              {
                                PushPtPair(Region._PtPairList, PtNbList[i], PtPtr, enSPHPtPair, -dxiac, -dyiac, driac2);
                              }
                          }
                        break;
                      }
						
                    case enNULLPt:
                      PushPtPair(Region._PtPairList, PtPtr, PtNbList[i], enSPHNULLPtPair, dxiac, dyiac, driac2);   
                      break;
                                                                               
                    case enDumPt://2014.09.07
                      PushPtPair(Region._PtPairList, PtPtr, PtNbList[i], enSPHDumPtPair, dxiac, dyiac, driac2);   
                      break;

                    case enBndPt:
                      PushPtPair(Region._PtPairList, PtPtr, PtNbList[i], enSPHBndPtPair, dxiac, dyiac, driac2);   
                      break;

                    case enEHDBndPt://2019.10.26
                      PushPtPair(Region._PtPairList, PtPtr, PtNbList[i], enSPHEHDBndPtPair, dxiac, dyiac, driac2);   
                      break;

                    case enEHDDumPt://2019.10.26
                      PushPtPair(Region._PtPairList, PtPtr, PtNbList[i], enSPHEHDDumPtPair, dxiac, dyiac, driac2);   
                      break;

                    default:
                      break;
                    }

                }
            }

          PtNbList.clear();
        }
    }

  //2019.10.26 search EHD Dummy neighours of EHD Boundary 
	for(size_t f=0;f!=Region._PtList.size(); ++f)
    {
      PtPtr=&Region._PtList[f];      

      if(PtPtr->_Type==enEHDBndPt)
        {
          _Box.GetNbl(PtPtr->_x,PtPtr->_y,PtPtr->_r,_PtNbList);

          for(size_t i=0;i<_PtNbList.size();i++)//for all the number of the _PtNbList
            {
              if(_PtNbList[i]->_Type==enEHDDumPt)
                {
                  dxiac=PtPtr->_x-_PtNbList[i]->_x;
                  dyiac=PtPtr->_y-_PtNbList[i]->_y;
                  driac2=dxiac*dxiac+dyiac*dyiac;
                  
                  PushPtPair(Region._PtPairList, PtPtr, _PtNbList[i], enEHDBndDumPtPair, dxiac, dyiac, driac2);  
                }
            }            
        }
    }
  //period boundary condtions && electric dummy particle neighbour searching 
  if(Region._ControlSPH._PerdBnd==1)
    {
      MAXX=Region._ControlSPH._PerdBndMaxX;
      MINX=Region._ControlSPH._PerdBndMinX;

      for(unsigned int f=0;f<Region._PtList.size();f++)
        {
          PtPtr=&Region._PtList[f];
          if(PtPtr->_Type==enEHDBndPt)
            {
              _Box.GetNbl2(PtPtr->_x, PtPtr->_y, PtPtr->_r,_PtNbList,Region);
				
              for(unsigned int i=0;i!=_PtNbList.size();++i)
                {
                  if(_PtNbList[i]->_Type==enEHDDumPt)
                    {
                      if(_PtNbList[i]->_x>PtPtr->_x)//主粒子在左侧边界附近
                        dxiac=PtPtr->_x-_PtNbList[i]->_x+(MAXX-MINX);
                      else//主粒子在右侧边界附近
                        dxiac=PtPtr->_x-_PtNbList[i]->_x-(MAXX-MINX);

                      dyiac=PtPtr->_y-_PtNbList[i]->_y;
                      dxiac2=dxiac*dxiac;
                      dyiac2=dyiac*dyiac;
                      driac2=dxiac2+dyiac2;
                      r=_PtNbList[i]->_r+PtPtr->_r;
                      r2=r*r;

                      PushPtPair(Region._PtPairList, PtPtr, _PtNbList[i], enEHDBndDumPtPair, dxiac, dyiac, driac2); 
                    }
                }
            }
        }
    }

  //to get the exact neighbour number for each SPH particle
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      Region._PtList[Region._CalList[i]]._NumNegbor--;
    }
  
  //表面张力粒子对搜索
  if(Region._ControlSPH._SPHST==1)		
    {
      for(size_t f=0;f<Region._CalList.size();f++)
        {
          ID=Region._CalList[f];

          PtPtr=&Region._PtList[ID];

          if(PtPtr->_Type==enSPHPt)
            {
              _Box.GetNbl(PtPtr->_x,PtPtr->_y,PtPtr->_rr,_PtNbList);

              for(size_t i=0;i<_PtNbList.size();i++)//for all the number of the _PtNbList
                {
                  dxiac=PtPtr->_x-_PtNbList[i]->_x;
                  dyiac=PtPtr->_y-_PtNbList[i]->_y;

                  dxiac2=dxiac*dxiac;
                  dyiac2=dyiac*dyiac;

                  driac2=dxiac2+dyiac2;

                  rr=_PtNbList[i]->_rr+Region._PtList[ID]._rr;
                  rr2=rr*rr;

                  //the second searching
                  if(4*driac2<=rr2)
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
					
                                Region._PtPairList2.push_back(PtPair);
                                break;
                              }

                            else
                              {
                                double PtNbr;
                                PtNbr=_PtNbList[i]->_rr*_PtNbList[i]->_rr;
                                if(driac2>PtNbr)
                                  {
                                    PtPair._Type=enSPHPtPair;
                                    PtPair._PtiPtr=&Region._PtList[ID];
                                    PtPair._PtjPtr=_PtNbList[i];
                                    PtPair._driac2=driac2;
						
                                    Region._PtPairList2.push_back(PtPair);
                                  }
                                break;
                              }
                          }
							
                        case enNULLPt: 
                          {
                            PtPair._Type=enSPHNULLPtPair;
                            PtPair._PtiPtr=&Region._PtList[ID];
                            PtPair._PtjPtr=_PtNbList[i];
                            PtPair._driac2=driac2;
						
                            Region._PtPairList2.push_back(PtPair);
                            break;
                          }

                        case enBndPt:
                          {
                            PtPair._Type=enSPHBndPtPair;
                            PtPair._PtiPtr=&Region._PtList[ID];
                            PtPair._PtjPtr=_PtNbList[i];
                            PtPair._driac2=driac2;
						
                            Region._PtPairList2.push_back(PtPair);
                            break;
                          }

                        default:
                          break;
                        }

                    }
                }
            }
        }
    }
}

void CNblSch::Clear(CRegion &Region)
{
	////delete P;
	//for(size_t j=0;j<Region._PtPairList.size();j++)
	//	if(Region._PtPairList[j]!=NULL)
	//	{
	//		delete Region._PtPairList[j];
	//		Region._PtPairList[j]=NULL;
	//	}
	Region._PtPairList.clear();

	//for(size_t j=0;j<Region._PtPairList2.size();j++)
	//if(Region._PtPairList2[j]!=NULL)
	//{
	//	delete Region._PtPairList2[j];
	//	Region._PtPairList2[j]=NULL;
	//}

	if(Region._ControlSPH._SPHST==1)
	{
		Region._PtPairList2.clear();
	}
}

void CNblSch::GetMshNbl( CRegion& Region )
{
	double dxiac;
	double dyiac;
	double dxiac2;
	double dyiac2;
	double driac2;
	double r,rr;
	double r2,rr2;
	unsigned int ID;
	CPtMshPair PtMshPair;

	_Box.UpdateBox(Region);//update the Box

	//周期性边界，目前只能用于X方向的周期性边界
	//周期性边界1.找出每一行的粒子坐标最大值和最小值
	if(Region._ControlSPH._PerdBnd==1)
	{
		vector<CBasePt *> BndPt;
		double MAXX,MINX;

		MAXX=Region._ControlSPH._PerdBndMaxX;
		MINX=Region._ControlSPH._PerdBndMinX;

		//周期性边界3.周期性边界的粒子对搜索
		for(unsigned int f=0;f<Region._MeshList.size();f++)
		{
			//if(Region._MeshList[f]._Type==enSPHPt)
			{
				_Box.GetNbl2(Region._MeshList[f]._x,Region._MeshList[f]._y,Region._MeshList[f]._r,_PtNbList,Region);

				for(unsigned int i=0;i!=_PtNbList.size();++i)
				{
					if(_PtNbList[i]->_x>Region._MeshList[f]._x)//主粒子在左侧边界附近
						dxiac=_PtNbList[i]->_x-Region._MeshList[f]._x-(MAXX-MINX);
					else//主粒子在右侧边界附近
						dxiac=_PtNbList[i]->_x-Region._MeshList[f]._x+(MAXX-MINX);

					dyiac=_PtNbList[i]->_y-Region._MeshList[f]._y;
					dxiac2=dxiac*dxiac;
					dyiac2=dyiac*dyiac;
					driac2=dxiac2+dyiac2;
					r=_PtNbList[i]->_r+Region._MeshList[f]._r;
					r2=r*r;

					//the second searching
					if(4*driac2<=r2)
					{
						PtMshPair._driac2=driac2;
						PtMshPair._MshjPtr=&Region._MeshList[f];
						PtMshPair._PtiPtr=_PtNbList[i];
						Region._PtMshPairList.push_back(PtMshPair);
					}
				}
			}
		}
	}


	//正常找粒子对部分
	//只用流体粒子作为主搜粒子
	for(size_t f=0;f<Region._MeshList.size();f++)
	{

		//if(Region._MeshList[ID]._Type==enSPHPt)
		{
			_Box.GetNbl(Region._MeshList[f]._x,Region._MeshList[f]._y,Region._MeshList[f]._r,_PtNbList);

			for(size_t i=0;i<_PtNbList.size();i++)//for all the number of the _PtNbList
			{
				dxiac=_PtNbList[i]->_x-Region._MeshList[f]._x;
				dyiac=_PtNbList[i]->_y-Region._MeshList[f]._y;

				dxiac2=dxiac*dxiac;
				dyiac2=dyiac*dyiac;

				driac2=dxiac2+dyiac2;
				r=_PtNbList[i]->_r+Region._MeshList[f]._r;
				r2=r*r;

				//the second searching
				if(4*driac2<=r2)
				{
					PtMshPair._MshjPtr=&Region._MeshList[f];
					PtMshPair._PtiPtr=_PtNbList[i];
					PtMshPair._driac2=driac2;

					Region._PtMshPairList.push_back(PtMshPair);
				}
			}
		}
	}

}

//used to push back particle pairs into PtPairList
void CNblSch::PushPtPair(vector<CPtPair> & PtPairList, CBasePt * PtiPtr, CBasePt * PtjPtr, enPTPAIRTYPE PtPairType,  double xij, double yij, double r2)
{
  CPtPair PtPair;
  PtPair._Type=PtPairType;
  PtPair._driac2=r2;
  PtPair._xij=xij;
  PtPair._yij=yij;
  PtPair._PtiPtr=PtiPtr;
  PtPair._PtjPtr=PtjPtr;
                              
  PtPairList.push_back(PtPair);

  //邻近粒子统计
  PtiPtr->_NumNegbor++;
  PtjPtr->_NumNegbor++;
}
