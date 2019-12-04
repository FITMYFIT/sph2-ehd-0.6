#include "Model.h"

CModel::CModel()
{
}

CModel::~CModel()
{
}

void CModel::EHDPlannar(CRegion &Region)
{
  unsigned int N=33;//the particle number in x and y directions
  //N=65;

  Region._ControlSPH._CellNumx=Region._ControlSPH._CellNumy=N;//the box number, specified here, the value in .k file will not be used

  double xlb=-0.5;//left bottom corner
  double ylb=-0.5;
  double xru=0.5;//right upper corner
  double yru=0.5;

  double psize=(xru-xlb)/(N);// particle size

  // Region._PtList.clear();
  Region._PtList.resize((N+1)*(N+1));

  CBasePt * BasePtPtr;
  for(unsigned int i=0;i!=N+1;++i)
    {
      for(unsigned int j=0;j!=N+1;++j)
        {
          // CBasePt * BasePtPtr=new(CBasePt);
          BasePtPtr=&Region._PtList[i*(N+1)+j];
          BasePtPtr->_ID=i*(N+1)+j+1;

          BasePtPtr->_x=xlb+j*psize;
          BasePtPtr->_y=ylb+i*psize;

          if(BasePtPtr->_y<0)
            {
              if(ISZERO((BasePtPtr->_y-ylb))) //bottom layer
                {
                  BasePtPtr->_PID=1;
                  BasePtPtr->_ePhi=1.0;
                  BasePtPtr->_Type=enEHDBndPt;
                }
              else
                {
                  BasePtPtr->_PID=2;
                  BasePtPtr->_ePhi=0.0;
                  BasePtPtr->_Type=enSPHPt;
                }
            }
          else
            {
              if(ISZERO((BasePtPtr->_y-yru))) //upper layer
                {
                  BasePtPtr->_PID=4;
                  BasePtPtr->_ePhi=0.0;
                  BasePtPtr->_Type=enEHDBndPt;
                }
              else
                {
                  BasePtPtr->_PID=3;
                  BasePtPtr->_ePhi=0.0;
                  BasePtPtr->_Type=enSPHPt;
                }
            }

          BasePtPtr->_Volume=psize*psize;

          BasePtPtr->_u=0.0;
          BasePtPtr->_v=0.0;

          BasePtPtr->_rho=1.0;
          BasePtPtr->_rho0=1.0;
          BasePtPtr->_m=BasePtPtr->_rho0*BasePtPtr->_Volume;
          BasePtPtr->_h=Region._PartList[BasePtPtr->_PID-1]._HdivDp*psize;
          BasePtPtr->_r=2.0*BasePtPtr->_h;

          BasePtPtr->_C0=Region._PartList[BasePtPtr->_PID-1]._C0;

          BasePtPtr->_Cs=Region._ControlSPH._Cs;

          BasePtPtr->_p=0.0;
          BasePtPtr->_Iflag=0;

          //ehd concerned variables
          BasePtPtr->_eEpsilon=0.0;
          BasePtPtr->_eEx=0.0;
          BasePtPtr->_eEy=0.0;
          BasePtPtr->_eRho=0.0;

          // Region._PtList.push_back(*BasePtPtr);
          // Region._PtList[i*N+j-1]=*BasePtPtr;
        }
    }

  // // //the planner model, add 2 layers upper and bottom
  unsigned int icount=0;
   unsigned int PtNum0=Region._PtList.size();
   CPart Part5,Part6;//the part on the bottom and uppercase
   Part5=Region._PartList[0];
   Part6=Region._PartList[3];   

   Part5._PID=5;
   Part5._PartType=enEHDDum;

   Part6._PID=6;
   Part6._PartType=enEHDDum;

   Region._PartList.push_back(Part5);
   Region._PartList.push_back(Part6);
   
   for (unsigned int i =0;i!=PtNum0;++i)
     {
       if(Region._PtList[i]._PID==1||Region._PtList[i]._PID==4)//the bottom layer
         {
           CBasePt BasePt2;
           CBasePt BasePt3;
           
           BasePt2=Region._PtList[i];
           BasePt3=Region._PtList[i];

           if(BasePt2._PID==1)//the bottom layer
             {
               BasePt2._y-=psize;
               BasePt3._y-=2*psize;
               BasePt2._PID=5;
               BasePt3._PID=5;
               BasePt2._ID=PtNum0+(++icount);
               BasePt3._ID=PtNum0+(++icount);
               
               BasePt2._Type=enEHDDumPt;
               BasePt3._Type=enEHDDumPt;
               
               Region._PtList.push_back(BasePt2);
               // Region._PtList.push_back(BasePt3);
             }
           if(BasePt2._PID==4)//the upper layer
             {
               BasePt2._y+=psize;
               BasePt3._y+=2*psize;
               BasePt2._PID=6;
               BasePt3._PID=6;
               BasePt2._ID=PtNum0+(++icount);
               BasePt3._ID=PtNum0+(++icount);

               BasePt2._Type=enEHDDumPt;
               BasePt3._Type=enEHDDumPt;
               
               Region._PtList.push_back(BasePt2);
               //Region._PtList.push_back(BasePt3);
             }
         }
     }

   //set ID of all the particles again
   icount=0;
   for(unsigned int i=0;i!=Region._PtList.size();++i)
     {
       icount++;
       Region._PtList[i]._ID=icount;
     }
}


void CModel::EHDBulkRelax(CRegion & Region)
{
  unsigned int N=64;//the particle number in x and y directions
  Region._ControlSPH._CellNumx=Region._ControlSPH._CellNumy=N;//the box number, specified here, the value in .k file will not be used

  double xlb=-0.5;//left bottom corner
  double ylb=-0.5;
  double xru=0.5;//right upper corner
  double yru=0.5;

  double psize=(xru-xlb)/(N);// particle size

  double xmin=xlb-psize;
  double ymin=ylb-psize;

  // Region._PtList.clear();
  Region._PtList.resize((N+3)*(N+3));

  double a=0.05;
  double r2;

  CBasePt * BasePtPtr;
  for(unsigned int i=0;i!=N+3;++i)
    {
      for(unsigned int j=0;j!=N+3;++j)
        {
          BasePtPtr=&Region._PtList[i*(N+3)+j];
          BasePtPtr->_ID=i*(N+3)+j+1;

          BasePtPtr->_x=xmin+j*psize;
          BasePtPtr->_y=ymin+i*psize;

          if(i==0||i==N+2||j==0||j==N+2)//the outer layer, end dummy pt
            {
              BasePtPtr->_PID=3;
              BasePtPtr->_Type=enEHDDumPt;
            }

          else if(i==1||i==N+1||j==1||j==N+1)
            {
              BasePtPtr->_PID=2;
              BasePtPtr->_Type=enEHDBndPt;
            }
          else
            {
              BasePtPtr->_PID=1;
              BasePtPtr->_Type=enSPHPt;
            }
       
          BasePtPtr->_Volume=psize*psize;

          BasePtPtr->_u=0.0;
          BasePtPtr->_v=0.0;

          BasePtPtr->_rho=1.0;
          BasePtPtr->_rho0=1.0;
          BasePtPtr->_m=BasePtPtr->_rho0*BasePtPtr->_Volume;
          BasePtPtr->_h=Region._PartList[BasePtPtr->_PID-1]._HdivDp*psize;
          BasePtPtr->_r=2.0*BasePtPtr->_h;

          BasePtPtr->_C0=Region._PartList[BasePtPtr->_PID-1]._C0;

          BasePtPtr->_Cs=Region._ControlSPH._Cs;

          BasePtPtr->_p=0.0;
          BasePtPtr->_Iflag=0;

          //ehd concerned variables
          BasePtPtr->_eEpsilon=0.0;
          BasePtPtr->_eEx=0.0;
          BasePtPtr->_eEy=0.0;
          BasePtPtr->_eRho=0.0;

          BasePtPtr->_ePhi=0.0;

          //for rhoe
          r2=pow(BasePtPtr->_x,2)+pow(BasePtPtr->_y,2);
          BasePtPtr->_eRho=pow(e,-r2/(2*a*a))/(a*sqrt(2*PI));

        }
    }

}


void CModel::EHDIsoCondCylinder (CRegion & Region)
{
  unsigned int N=64;//the particle number in x and y directions
  N=128;
  Region._ControlSPH._CellNumx=Region._ControlSPH._CellNumy=N;//the box number, specified here, the value in .k file will not be used

  double xlb=-0.5;//left bottom corner
  double ylb=-0.5;
  double xru=0.5;//right upper corner
  double yru=0.5;

  double psize=(xru-xlb)/(N);// particle size

  double xmin=xlb-psize;
  double ymin=ylb-psize;

  // Region._PtList.clear();
  Region._PtList.resize((N+3)*(N+3));

  double a=0.05;
  double r2;
  double R=0.05;
  CBasePt * BasePtPtr;
  for(unsigned int i=0;i!=N+3;++i)
    {
      for(unsigned int j=0;j!=N+3;++j)
        {
          BasePtPtr=&Region._PtList[i*(N+3)+j];
          BasePtPtr->_ID=i*(N+3)+j+1;

          BasePtPtr->_x=xmin+j*psize;
          BasePtPtr->_y=ymin+i*psize;

          if(i==0||i==N+2||j==0||j==N+2)//the outer layer, end dummy pt
            {
              BasePtPtr->_PID=4;
              BasePtPtr->_Type=enEHDDumPt;
              BasePtPtr->_eRho=0.5;
            }

          else if(i==1||i==N+1||j==1||j==N+1)
            {
              BasePtPtr->_PID=3;
              BasePtPtr->_Type=enEHDBndPt;
              BasePtPtr->_eRho=0.5;
            }
          else
            {
              BasePtPtr->_Type=enSPHPt;
              if(pow(BasePtPtr->_x,2)+pow(BasePtPtr->_y,2)<=R*R)
                {
                  BasePtPtr->_PID=1;
                  BasePtPtr->_eRho=0.5;
                }
              else
                {
                  BasePtPtr->_PID=2;
                  BasePtPtr->_eRho=0;
                }            
            }
       
          BasePtPtr->_Volume=psize*psize;

          BasePtPtr->_u=0.0;
          BasePtPtr->_v=0.0;

          BasePtPtr->_rho=1.0;
          BasePtPtr->_rho0=1.0;
          BasePtPtr->_m=BasePtPtr->_rho0*BasePtPtr->_Volume;
          BasePtPtr->_h=Region._PartList[BasePtPtr->_PID-1]._HdivDp*psize;
          BasePtPtr->_r=2.0*BasePtPtr->_h;

          BasePtPtr->_C0=Region._PartList[BasePtPtr->_PID-1]._C0;

          BasePtPtr->_Cs=Region._ControlSPH._Cs;

          BasePtPtr->_p=0.0;
          BasePtPtr->_Iflag=0;

          //ehd concerned variables
          BasePtPtr->_eEpsilon=0.0;
          BasePtPtr->_eEx=0.0;
          BasePtPtr->_eEy=0.0;
          //BasePtPtr->_eRho=0.0;

          BasePtPtr->_ePhi=0.0;

          }
    }

}


void CModel::EHDDrop (CRegion & Region)
{
  unsigned int N=64;//the particle number in x and y directions
  // N=128;
 N=256;
 // N=512;
 unsigned int Nx=N;
 unsigned int Ny=N;
 Region._ControlSPH._CellNumx=Region._ControlSPH._CellNumy=N;//the box number, specified here, the value in .k file will not be used

 double xlb=-1;//left bottom corner
 double ylb=-1;
 double xru=1;//right upper corner
 double yru=1;

 double psize=(xru-xlb)/(N);// particle size

 double xmin=xlb-psize;
 double ymin=ylb-psize;

 // Region._PtList.clear();
 Region._PtList.clear();
 Region._PtList.resize((Ny+3)*(Nx+3));

 double r2;
 double R0=0.1;
 CBasePt * BasePtPtr;
  
 for(unsigned int i=0;i!=Ny+3;++i)
     {
       for(unsigned int j=0;j!=Nx+3;++j)
         {
           BasePtPtr=&Region._PtList[i*(Nx+3)+j];
           BasePtPtr->_ID=i*(Nx+3)+j+1;

           BasePtPtr->_x=xmin+j*psize;
           BasePtPtr->_y=ymin+i*psize;

           if(i==0||i==Ny+2||j==0||j==Nx+2)//the outer layer, end dummy pt
             //if(i==0||i==Ny+2)//the outer layer, ehd dummy pt
             {
               BasePtPtr->_PID=4;
               BasePtPtr->_Type=enEHDDumPt;
               BasePtPtr->_eRho=0.0;
             }

           else if(i==1||i==Ny+1||j==1||j==Nx+1)
             // else if(i==1||i==Ny+1)
             {
               BasePtPtr->_PID=3;
               BasePtPtr->_Type=enEHDBndPt;
               BasePtPtr->_eRho=0.0;
             }
           else
             {
               BasePtPtr->_Type=enSPHPt;
               if(pow(BasePtPtr->_x,2)+pow(BasePtPtr->_y,2)<=R0*R0)
                 //if(fabs(BasePtPtr->_x)<=R0&&fabs(BasePtPtr->_y)<=R0)
                 {
                   BasePtPtr->_PID=1;
                   BasePtPtr->_eRho=0.0;
                 }
               else
                 {
                   BasePtPtr->_PID=2;
                   BasePtPtr->_eRho=0;
                 }            
             }
       
           BasePtPtr->_Volume=psize*psize;

           BasePtPtr->_u=0.0;
           BasePtPtr->_v=0.0;

           BasePtPtr->_rho=1.0;
           BasePtPtr->_rho0=1.0;
           BasePtPtr->_m=BasePtPtr->_rho0*BasePtPtr->_Volume;
           BasePtPtr->_h=Region._PartList[BasePtPtr->_PID-1]._HdivDp*psize;
           BasePtPtr->_r=2.0*BasePtPtr->_h;

           BasePtPtr->_C0=Region._PartList[BasePtPtr->_PID-1]._C0;

           BasePtPtr->_Cs=Region._ControlSPH._Cs;

           BasePtPtr->_p=0.0;
           BasePtPtr->_Iflag=0;

           //ehd concerned variables
           BasePtPtr->_eEpsilon=0.0;
           BasePtPtr->_eEx=0.0;
           BasePtPtr->_eEy=0.0;
           //BasePtPtr->_eRho=0.0;

           if(Region._ControlSPH._SPHST==1)
             {
               BasePtPtr->_H=Region._PartList[BasePtPtr->_PID-1]._STHvsh*BasePtPtr->_h;
               BasePtPtr->_rr=2*BasePtPtr->_H;
             }


           double Einf=1.0;//4.3 step 1 ,test for static case
           if(BasePtPtr->_PID!=1)
             BasePtPtr->_ePhi=Einf*BasePtPtr->_x;

         }
     }


}

void CModel::EHDDrop2 (CRegion & Region)//input model from truegrid model
{
  unsigned int N=64;//the particle number in x and y directions
  N=128;
  // N=256;
  // N=512;
  unsigned int Nx=N;
  unsigned int Ny=N;
  Region._ControlSPH._CellNumx=Region._ControlSPH._CellNumy=N;//the box number, specified here, the value in .k file will not be used

  double xlb=-1;//left bottom corner
  double ylb=-1;
  double xru=1;//right upper corner
  double yru=1;

  double psize=(xru-xlb)/(N);// particle size

  double xmin=xlb-psize;
  double ymin=ylb-psize;

  // Region._PtList.clear();
 
  double r2;
  double R0=0.1;
  CBasePt * BasePtPtr;
  
  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      BasePtPtr=&Region._PtList[i];
      
      //ehd concerned variables
      BasePtPtr->_eEpsilon=0.0;
      BasePtPtr->_eEx=0.0;
      BasePtPtr->_eEy=0.0;
      //BasePtPtr->_eRho=0.0;

      double Einf=1.0;//4.3 step 1 ,test for static case
      if(BasePtPtr->_PID!=1)
        BasePtPtr->_ePhi=Einf*BasePtPtr->_x;

  
    }


}
