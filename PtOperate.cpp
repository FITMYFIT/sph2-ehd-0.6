#include "PtOperate.h"


double PtDistance(CBasePt * PtiPtr, CBasePt * PtjPtr)//the distance of 2 particles
{
  double x1,y1;
  double x2,y2;
  double dis;

  x1=PtiPtr->_x;
  y1=PtiPtr->_y;

  x2=PtjPtr->_x;
  y2=PtjPtr->_y;

  dis=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

  return dis;
}

double PtDistance2(CBasePt * PtiPtr, CBasePt * PtjPtr)//the distance of 2 particles
{
  double x1,y1;
  double x2,y2;
  double dis;

  x1=PtiPtr->_x;
  y1=PtiPtr->_y;

  x2=PtjPtr->_x;
  y2=PtjPtr->_y;

  dis=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);

  return dis;
}
