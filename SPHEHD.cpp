#include "SPHEHD.h"


CSPHEHD::CSPHEHD()
{
}

CSPHEHD::~CSPHEHD()
{
}

void CSPHEHD::Solve(CRegion & Region, unsigned int Timesteps)
{
  CSPHPt * PtiPtr,*PtjPtr;
  CBasePt * BasePtPtr;
  CKnl * KnlPtr;
  double Exi,Exj;
  double Eyi,Eyj;
  double epsiloni, epsilonj;//epsilon of part
  double kappai, kappaj;//kappa of part
  double norm;
  double normx,normy;
  double normx1,normy1;
  unsigned int IDi,IDj,ID2;

  // CKFile kfiletemp;
  // kfiletemp.outTecplotEHDPLANNER(Region,Timesteps,"xxehd-case1");
  // kfiletemp.OutPtNeighbour(Region, &Region._PtList[10], Timesteps, "neighbour10");

  //interpolate Ex,Ey, Phie of ehd boundary particle from fluid particles
  InterpEHDBndProperty(Region, Timesteps);

  //intepolate phi of ehd dummy particles from fluid and ehd boundary particles
  InterpEHDDumPropterty(Region, Timesteps);
  
  // //1.update eepsilon and ekappa
  //Lopez 2011 use an interpolation scheme to calculate transition value between the two fluid with VOF colour
  //here we try an direct interpolation using SPH, just like the interpolation of colour
  for (unsigned int i=0;i!=Region._PtPairList.size();i++)
    {
      // if(Region._PtPairList[i]._Type!=enSphbndptpair&&Region._PtPairList[i]._Type!=enSPHDumPtPair)
      {
        PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
        PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
        KnlPtr=&Region._KnlList[i];

        epsiloni=Region._PartList[PtiPtr->_PID-1]._eEpsilon;
        epsilonj=Region._PartList[PtjPtr->_PID-1]._eEpsilon;

        kappai=Region._PartList[PtiPtr->_PID-1]._eKappa;
        kappaj=Region._PartList[PtjPtr->_PID-1]._eKappa;

        if(PtiPtr->_Type==enSPHPt)
          {
            PtiPtr->_eEpsilon+=PtjPtr->_mrho*KnlPtr->_W*epsilonj;

            PtiPtr->_eKappa+=PtjPtr->_mrho*KnlPtr->_W*kappaj;
          }
        
        if(PtjPtr->_Type==enSPHPt)         
          {
            if(PtiPtr!=PtjPtr)
              {
                PtjPtr->_eEpsilon+=PtiPtr->_mrho*KnlPtr->_W*epsiloni;

                PtjPtr->_eKappa+=PtiPtr->_mrho*KnlPtr->_W*kappai;
              }
          }

        //for ehd dummy or ehd boundary particles
        else
          {
            PtjPtr->_eEpsilon=epsilonj;
            PtjPtr->_eKappa=kappaj;
          }
      }
    }
  
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //temp, specify epsilon and kappa from the part, without transient region between two fluid
  // for(unsigned int i=0;i!=Region._PtList.size();++i)
  // {
  //   BasePtPtr=&Region._PtList[i];
  //   BasePtPtr->_eEpsilon=Region._PartList[BasePtPtr->_PID-1]._eEpsilon;
  // }

  
  //2.solve charge density conservation equ Lopez 2011 Equ(21)
  //2.1 solve deRho
  for (unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      // if(Region._PtPairList[i]._Type==enSPHPtPair)
      {
        PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
        PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
        KnlPtr=&Region._KnlList[i];

        if(Region._ControlSPH._RunMod==1)//Explicit
          {

            if(PtiPtr->_Type==enSPHPt)
              {
                PtiPtr->_deRho-=PtjPtr->_mrho*((PtiPtr->_eKappa*PtjPtr->_eEx+PtjPtr->_eKappa*PtiPtr->_eEx)*KnlPtr->_Wx
                                               +(PtiPtr->_eKappa*PtjPtr->_eEy+PtjPtr->_eKappa*PtiPtr->_eEy)*KnlPtr->_Wy);
                PtiPtr->_deRho-=PtjPtr->_mrho*((PtiPtr->_eRho*PtjPtr->_u+PtjPtr->_eRho*PtiPtr->_u)*KnlPtr->_Wx
                                               +(PtiPtr->_eRho*PtjPtr->_v+PtjPtr->_eRho*PtiPtr->_v)*KnlPtr->_Wy);

                //PtiPtr->_deRho-=PtjPtr->_mrho*PtjPtr->_eKappa*(PtjPtr->_eEx*KnlPtr->_Wx+PtjPtr->_eEy*KnlPtr->_Wy);
              }
            if((PtiPtr!=PtjPtr)&&PtjPtr->_Type==enSPHPt)
              {
                PtjPtr->_deRho+=PtiPtr->_mrho*((PtjPtr->_eKappa*PtiPtr->_eEx+PtiPtr->_eKappa*PtjPtr->_eEx)*KnlPtr->_Wx
                                               +(PtjPtr->_eKappa*PtiPtr->_eEy+PtiPtr->_eKappa*PtjPtr->_eEy)*KnlPtr->_Wy);
                PtjPtr->_deRho+=PtiPtr->_mrho*((PtjPtr->_eRho*PtiPtr->_u+PtiPtr->_eRho*PtjPtr->_u)*KnlPtr->_Wx
                                               +(PtjPtr->_eRho*PtiPtr->_v+PtiPtr->_eRho*PtjPtr->_v)*KnlPtr->_Wy);

                // PtjPtr->_deRho+=PtiPtr->_mrho*PtjPtr->_eKappa*(PtiPtr->_eEx*KnlPtr->_Wx+PtiPtr->_eEy*KnlPtr->_Wy);
              }
          }

        else//Implicit
          {
            IDi=PtiPtr->_ID2;
          }
      }
    }

  //2.2 update volume charge density eRho, rhoe n+1
  double DeltaT;
  DeltaT=Region._ControlSPH._DeltaT;
  for(unsigned int i=0; i<Region._PtList.size();i++)
    {
      BasePtPtr=&Region._PtList[i];
      if(BasePtPtr->_Type==enSPHPt)
        {
          BasePtPtr->_eRho+=DeltaT*BasePtPtr->_deRho;
        }
    }

  //3. calculate electric field E
  //3.1 calculate φ, Lopez2011 Equ.(22),Poisson equ.  
  //3.1.1 define and initialize the matrix A and right hand side b, here I use LASpack
  //---------------------------------------------------------------------------------------------------------------------
  unsigned int PtNum=Region._CalList.size();
  //use laspack to solve linear equ.
  QMatrix LasA;//laspack amtrix A
  Vector Lasx,Lasb;//laspack vector, not std vector
  size_t Dim=PtNum;

  char namea[]="A";
  char namex[]="x";
  char nameb[]="b";
  Q_Constr(&LasA, namea, Dim , False , Rowws, Normal , True );
  V_Constr(&Lasx, namex, Dim, Normal, True);
  V_Constr(&Lasb, nameb, Dim, Normal, True);

  //set the length of LasA, after q_setlen, element will be set to 0
  //for i, some neighbor will not involved in calculation (if neighbor j is dumpt or not in calculate range)
  //so the length of row may larger than the element number, but it does not influence the result of Ax=b
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];

      Q_SetLen(&LasA, BasePtPtr->_ID2, BasePtPtr->_NumNegbor);

      Q_SetEntry(&LasA, BasePtPtr->_ID2, 0, BasePtPtr->_ID2, 0.0);//set the first one of row as i-i element
    }
  //set all components of lasvector 0
  V_SetAllCmp(&Lasb, 0.0);

  //3.1.1 done
  //---------------------------------------------------------------------------------------------------------------------

  //3.1.2 contruct and solove Ax=b
  //\[ - {\left( {{\rho _e}} \right)_a} = \sum\limits_b {\frac{{{m_b}}}{{{\rho _b}}}} ({\varepsilon _a} + {\varepsilon _b})({\phi _a} - {\phi _b})\frac{1}{{{r_{ab}}}}\frac{{{\partial _a}{W_{ab}}}}{{\partial {r_{ab}}}}\]
  //the source term, b
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
      V_AddCmp(&Lasb, BasePtPtr->_ID2, -BasePtPtr->_eRho);
    }

   vector<unsigned int> icount;//element number of each row
   icount.resize(PtNum,1);//there is one element now, i-i, the first one of each row

   for(unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          //as (phi i-phi j ) is calculated,  so ptiptr=ptjptr makes no sense
          if(PtiPtr!=PtjPtr)
            {
              IDi=PtiPtr->_ID2;
              IDj=PtjPtr->_ID2;

              norm=(PtiPtr->_eEpsilon+PtjPtr->_eEpsilon)*KnlPtr->_Ww;
              // norm=4*(PtiPtr->_eEpsilon*PtjPtr->_eEpsilon)/(PtiPtr->_eEpsilon+PtjPtr->_eEpsilon)*KnlPtr->_Ww;//Monaghan 2005 Equ(7.4)

              Q_AddVal(&LasA, IDi, 0, PtjPtr->_mrho*norm);//the first one in row is i i

              if(PtjPtr->_Type==enSPHPt)
                {
                  icount[IDi-1]++;

                  Q_SetEntry(&LasA, IDi, icount[IDi-1]-1, IDj, -PtjPtr->_mrho*norm);

                  icount[IDj-1]++;

                  Q_AddVal(&LasA, IDj, 0, PtiPtr->_mrho*norm);//the first one in row j-j
                  Q_SetEntry(&LasA, IDj, icount[IDj-1]-1, IDi, -PtiPtr->_mrho*norm);

                }

              //ptjptr is a ehd dummy or ehd boundary particle
              else
                {
                  V_AddCmp(&Lasb, IDi, PtjPtr->_mrho*norm*PtjPtr->_ePhi);
                }
            }
        }
    }
  icount.clear();

  SetRTCAccuracy(1e-8);
  V_SetAllCmp(&Lasx, 0.0);
  // CGIter(&LasA, &Lasx, &Lasb, 100, SSORPrecond,1.2);
  // CGIter(&LasA, &Lasx, &Lasb, 20, SSORPrecond,1.2);
  // BiCGIter(&LasA, &Lasx, &Lasb, 100, SSORPrecond,1.2);
  // BiCGIter(&LasA, &Lasx, &Lasb, 20, SSORPrecond,1.2);
  // BiCGSTABIter(&LasA , &Lasx, &Lasb, 100, SSORPrecond,1.2);
  BiCGSTABIter(&LasA, &Lasx, &Lasb, 100, SSORPrecond,1.2);//the 4th parameter 100, is the maxiter steps
  // outputV( Lasx,  "xxbicgs-20");
  //output2(LasA, Lasx, Lasb, "xxbicgs-ax=b");

  //set φ of particle from b
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];

      BasePtPtr->_ePhi=V_GetCmp(&Lasx, BasePtPtr->_ID2);
    }

  Q_Destr(&LasA);
  V_Destr(&Lasx);
  V_Destr(&Lasb);

  //3.1.2 done
  //3.1 done
  //---------------------------------------------------------------------------------------------------------------------
  //3.2 calculate E using CSPM, as the dirct scheme produce unacceptable errors
  //3.2.1 cspm coefficient matrix
  if(Region._ControlSPH._CSPMIflag2==0)
    {
      _CSPMEHD.GetCSPMGradCorctCoef2(Region);//dummy particle involved
    }

  // //3.2.2 cspm source term for E=grad(phi)
  vector<double> bx;//CSPM修正的源项
	vector<double> by;
	bx.resize(Region._CalList.size(),0.0);
	by.resize(Region._CalList.size(),0.0);
  for(unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          IDi=PtiPtr->_ID2;
          IDj=PtjPtr->_ID2;

          if(PtiPtr!=PtjPtr)
            {
              bx[IDi-1]+=PtjPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wx;
              by[IDi-1]+=PtjPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wy;

              //particle j
              if(PtjPtr->_Type==enSPHPt)
                {
                  bx[IDj-1]+=PtiPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wx;
                  by[IDj-1]+=PtiPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wy;
                }
            }
        }
    }

  //calculate CSPM corrected E=grad(phi)
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
      ID2=BasePtPtr->_ID2;

      BasePtPtr->_eEx=-(BasePtPtr->_CSPMAxx*bx[ID2-1]+BasePtPtr->_CSPMAyx*by[ID2-1]);
      BasePtPtr->_eEy=-(BasePtPtr->_CSPMAxy*bx[ID2-1]+BasePtPtr->_CSPMAyy*by[ID2-1]);
    }

  bx.clear();
  by.clear();


  //3.3 calculate Fex,Fey:the accelerate caused by ehd force
  //fomulation details are in Richard's STUDY DIRARY "2019年9月13日星期五" part 4 (similar formulation with Lopez 2011 Equ32)
  //use CSPM
  //---------------------------------------------------------------------------------------------------------------------
  //3.2 calculate electric force Fe using CSPM
  //3.2.1 cspm coefficient matrix
  if(Region._ControlSPH._CSPMIflag2!=1)
    {
      _CSPMEHD.GetCSPMGradCorctCoef(Region);//dummy particle involved
    }

  //cspm source term
	bx.resize(Region._CalList.size(),0.0);
	by.resize(Region._CalList.size(),0.0);
  for(unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         // ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         // ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair
         )
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          IDi=PtiPtr->_ID2;
          IDj=PtjPtr->_ID2;
          
          epsiloni=PtiPtr->_eEpsilon;
          epsilonj=PtjPtr->_eEpsilon;

          Exi=PtiPtr->_eEx;
          Eyi=PtiPtr->_eEy;
          Exj=PtjPtr->_eEx;
          Eyj=PtjPtr->_eEy;

          normx=epsiloni*(Exi*(Exi*KnlPtr->_Wx+Eyi*KnlPtr->_Wy)-0.5*(Exi*Exi+Eyi*Eyi)*KnlPtr->_Wx)
            -epsilonj*(Exj*(Exj*KnlPtr->_Wx+Eyj*KnlPtr->_Wy)-0.5*(Exj*Exj+Eyj*Eyj)*KnlPtr->_Wx);
          normy=epsiloni*(Eyi*(Exi*KnlPtr->_Wx+Eyi*KnlPtr->_Wy)-0.5*(Exi*Exi+Eyi*Eyi)*KnlPtr->_Wy)
            -epsilonj*(Eyj*(Exj*KnlPtr->_Wx+Eyj*KnlPtr->_Wy)-0.5*(Exj*Exj+Eyj*Eyj)*KnlPtr->_Wy);

          if(PtiPtr!=PtjPtr)
            {
              bx[IDi-1]+=PtjPtr->_mrho*normx;
              by[IDi-1]+=PtjPtr->_mrho*normy;

              //particle j
              if(PtjPtr->_Type==enSPHPt)
                {
                  bx[IDj-1]+=PtiPtr->_mrho*normx;
                  by[IDj-1]+=PtiPtr->_mrho*normy;
                }
            }
        }
    }

  //calculate CSPM corrected electric force
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
      ID2=BasePtPtr->_ID2;

      BasePtPtr->_Fex=BasePtPtr->_CSPMAxx*bx[ID2-1]+BasePtPtr->_CSPMAyx*by[ID2-1];
      BasePtPtr->_Fey=BasePtPtr->_CSPMAxy*bx[ID2-1]+BasePtPtr->_CSPMAyy*by[ID2-1];
    }

  bx.clear();
  by.clear();

  CKFile kfiletemp;
  kfiletemp.outTecplotEHDPLANNER(Region,Timesteps,"xxehd-case2");
  // kfiletemp.OutPtNeighbour(Region, &Region._PtList[1290], Timesteps, "1291neighbor");

  std::cout<<"--------------------------------"<<endl;
}

// the scheme used in Basilisk,http://basilisk.dalembert.upmc.fr/src/ehd/implicit.h
void CSPHEHD::Solve2(CRegion & Region, unsigned int Timesteps)
{
  CSPHPt * PtiPtr,*PtjPtr;
  CBasePt * BasePtPtr;
  CKnl * KnlPtr;
  double epsiloni, epsilonj;//epsilon of part
  double kappai, kappaj;//kappa of part
  double norm;
  double normx,normy;

  double normx1,normy1;
  double Ex,Ey;
  double Exi,Exj,Eyi,Eyj;
  unsigned int IDi,IDj;  
  unsigned int ID,ID2;
  //1.update eepsilon and ekappa
  //Lopez 2011 use an interpolation scheme to calculate transition value between the two fluid with VOF colour
  //here we try an direct interpolation using SPH, just like the interpolation of colour
  for (unsigned int i=0;i!=Region._PtPairList.size();i++)
    {
      // if(Region._PtPairList[i]._Type!=enSphbndptpair&&Region._PtPairList[i]._Type!=enSPHDumPtPair)
      {
        PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
        PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
        KnlPtr=&Region._KnlList[i];

        epsiloni=Region._PartList[PtiPtr->_PID-1]._eEpsilon;
        epsilonj=Region._PartList[PtjPtr->_PID-1]._eEpsilon;

        kappai=Region._PartList[PtiPtr->_PID-1]._eKappa;
        kappaj=Region._PartList[PtjPtr->_PID-1]._eKappa;

        if(PtiPtr->_Type==enSPHPt)
          {
            PtiPtr->_eEpsilon+=PtjPtr->_mrho*KnlPtr->_W*epsilonj;

            PtiPtr->_eKappa+=PtjPtr->_mrho*KnlPtr->_W*kappaj;
          }
        
        if(PtjPtr->_Type==enSPHPt)         
          {
            if(PtiPtr!=PtjPtr)
              {
                PtjPtr->_eEpsilon+=PtiPtr->_mrho*KnlPtr->_W*epsiloni;

                PtjPtr->_eKappa+=PtiPtr->_mrho*KnlPtr->_W*kappai;
              }
          }

        //for ehd dummy or ehd boundary particles
        else
          {
            PtjPtr->_eEpsilon=epsilonj;
            PtjPtr->_eKappa=kappaj;
          }
      }
    }

  //temp, for test, do not interpolate kappa and epsilon
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //temp, specify epsilon and kappa from the part, without transient region between two fluid
  // for(unsigned int i=0;i!=Region._PtList.size();++i)
  //   {
  //     BasePtPtr=&Region._PtList[i];
  //     BasePtPtr->_eEpsilon=Region._PartList[BasePtPtr->_PID-1]._eEpsilon;
  //     BasePtPtr->_eKappa=Region._PartList[BasePtPtr->_PID-1]._eKappa;
  //   }
  
  //interpolate Ex,Ey, Phie of ehd boundary particle from fluid particles
  InterpEHDBndProperty(Region, Timesteps);

  //intepolate phi of ehd dummy particles from fluid and ehd boundary particles
  InterpEHDDumPropterty(Region, Timesteps);

  //2.calulate phi n+1, construct A(phi n+1)=b, nabla((K*deltaT+epsilon)*nabla(phi n+1))=-rhoe n
  //---------------------------------------------------------------------------------------------------------------------
  //use laspack to solve linear equ.
  unsigned int PtNum=Region._CalList.size();//the number used to initializae the size of linear equation
  QMatrix LasA;//laspack amtrix A
  Vector Lasx,Lasb;//laspack vector, not std vector
  size_t Dim=PtNum;

  char namea[]="A";
  char namex[]="x";
  char nameb[]="b";
  Q_Constr(&LasA, namea, Dim , False , Rowws, Normal , True );
  V_Constr(&Lasx, namex, Dim, Normal, True);
  V_Constr(&Lasb, nameb, Dim, Normal, True);

  //set the length of LasA, after q_setlen, element will be set to 0
  //for i, some neighbor will not involved in calculation (if neighbor j is dumpt or not in calculate range)
  //so the length of row may larger than the element number, but it does not influence the result of Ax=b
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];

      Q_SetLen(&LasA, BasePtPtr->_ID2, BasePtPtr->_NumNegbor);

      Q_SetEntry(&LasA, BasePtPtr->_ID2, 0, BasePtPtr->_ID2, 0.0);//set the first one of row as i-i element
    }
  //set all components of lasvector 0
  V_SetAllCmp(&Lasb, 0.0);

  //initilize laspack matrix A and vector b, done
  //---------------------------------------------------------------------------------------------------------------------

  //3.1.2 contruct and solove Ax=b
  // construct A(phi n+1)=b, nabla((K*deltaT+epsilon)*nabla(phi n+1))=-rhoe n
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
      V_AddCmp(&Lasb, BasePtPtr->_ID2, -BasePtPtr->_eRho);
    }

  vector<unsigned int> icount;//element number of each row
  icount.resize(PtNum,1);//there is one element now, i-i, the first one of each row
  double DeltaT;
  DeltaT=Region._ControlSPH._DeltaT;
  for(unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          //as (phi i-phi j ) is calculated,  so ptiptr=ptjptr makes no sense
          if(PtiPtr!=PtjPtr)
            {
              IDi=PtiPtr->_ID2;
              IDj=PtjPtr->_ID2;

              // epsiloni=Region._PartList[PtiPtr->_PID-1]._eEpsilon;//epsilon not smoothed
              // epsilonj=Region._PartList[PtjPtr->_PID-1]._eEpsilon;

              // kappai=Region._PartList[PtiPtr->_PID-1]._eKappa;//kappa not smoothed
              // kappaj=Region._PartList[PtjPtr->_PID-1]._eKappa;

              epsiloni=PtiPtr->_eEpsilon;
              epsilonj=PtjPtr->_eEpsilon;

              kappai=PtiPtr->_eKappa;
              kappaj=PtjPtr->_eKappa;

              norm=(kappai*DeltaT+epsiloni+kappaj*DeltaT+epsilonj)*KnlPtr->_Ww;
              //norm=2*(kappai*DeltaT+epsiloni)*KnlPtr->_Ww;//test for jump between two parts
              
              Q_AddVal(&LasA, IDi, 0, PtjPtr->_mrho*norm);//the first one in row is i i

              if(PtjPtr->_Type==enSPHPt)
                {
                  icount[IDi-1]++;

                  Q_SetEntry(&LasA, IDi, icount[IDi-1]-1, IDj, -PtjPtr->_mrho*norm);

                  icount[IDj-1]++;

                  //norm=2*(kappaj*DeltaT+epsilonj)*KnlPtr->_Ww;//test for jump between two parts
                  
                  Q_AddVal(&LasA, IDj, 0, PtiPtr->_mrho*norm);//the first one in row j-j
                  Q_SetEntry(&LasA, IDj, icount[IDj-1]-1, IDi, -PtiPtr->_mrho*norm);

                }

              //ptjptr is a ehd dummy or ehd boundary particle
              else
                {
                  V_AddCmp(&Lasb, IDi, PtjPtr->_mrho*norm*PtjPtr->_ePhi);
                }
            }
        }
    }
  icount.clear();

  SetRTCAccuracy(1e-8);
  //set the initial value using the previous timestep value
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
      V_SetCmp(&Lasx, BasePtPtr->_ID2, BasePtPtr->_ePhi);
    }
  // V_SetAllCmp(&Lasx, 0.0);//initial value 0

  BiCGSTABIter(&LasA, &Lasx, &Lasb, 100, SSORPrecond,1.2);//the 4th parameter 100, is the maxiter steps
  // outputV( Lasx,  "xxbicgs-20");
  //output2(LasA, Lasx, Lasb, "xxbicgs-ax=b");

  //set φ of particle from b
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];

      BasePtPtr->_ePhi=V_GetCmp(&Lasx, BasePtPtr->_ID2);
    }

  Q_Destr(&LasA);
  V_Destr(&Lasx);
  V_Destr(&Lasb);

  //rhoe n+1=-nabla*(epsilon*nabla(phi n+1))
  //http://basilisk.dalembert.upmc.fr/src/ehd/implicit.h
  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      Region._PtList[i]._eRho=0.0;
    }
  for(unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          //rhoei=-sum(mj/rhoj*(epsiloni+epsilonj)*(phii-phij)*Ww)
          if(PtiPtr!=PtjPtr)
            {
              // epsiloni=Region._PartList[PtiPtr->_PID-1]._eEpsilon;//epsilon not smoothed
              // epsilonj=Region._PartList[PtjPtr->_PID-1]._eEpsilon;

              epsiloni=PtiPtr->_eEpsilon;
              epsilonj=PtjPtr->_eEpsilon;

              norm=(epsiloni+epsilonj)*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Ww;
              // norm=2*epsiloni*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Ww;//test for jump between two parts
              
              PtiPtr->_eRho+=-PtjPtr->_mrho*norm;

              if(PtjPtr->_Type==enSPHPt)
                {
                  // norm=2*epsilonj*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Ww;//test for jump between two parts
                  
                  PtjPtr->_eRho+=PtiPtr->_mrho*norm;
                }            
            }
        }
    }
  
 
  //calculate E n+1=-grad(phi n+1)
  //---------------------------------------------------------------------------------------------------------------------
  //3.2 calculate E using CSPM, as the dirct scheme produce unacceptable errors
  //3.2.1 cspm coefficient matrix
  if(Region._ControlSPH._CSPMIflag2==0)
    {
      _CSPMEHD.GetCSPMGradCorctCoef2(Region);//dummy particle involved
    }

  // //3.2.2 cspm source term for E=grad(phi)
  vector<double> bx;//CSPM修正的源项
  vector<double> by;
  bx.resize(Region._CalList.size(),0.0);
  by.resize(Region._CalList.size(),0.0);
  for(unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          IDi=PtiPtr->_ID2;
          IDj=PtjPtr->_ID2;

          if(PtiPtr!=PtjPtr)
            {
              bx[IDi-1]+=PtjPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wx;
              by[IDi-1]+=PtjPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wy;

              //particle j
              if(PtjPtr->_Type==enSPHPt)
                {
                  bx[IDj-1]+=PtiPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wx;
                  by[IDj-1]+=PtiPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wy;
                }
            }
        }
    }

  //calculate CSPM corrected E=grad(phi)
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
      ID2=BasePtPtr->_ID2;

      BasePtPtr->_eEx=-(BasePtPtr->_CSPMAxx*bx[ID2-1]+BasePtPtr->_CSPMAyx*by[ID2-1]);
      BasePtPtr->_eEy=-(BasePtPtr->_CSPMAxy*bx[ID2-1]+BasePtPtr->_CSPMAyy*by[ID2-1]);
    }

  bx.clear();
  by.clear();
  //---------------------------------------------------------------------------

  //3.3 calculate the accelerate caused by ehd force
  //use equ 11, since the derivative of E is unaccruate and formulation of Lopez 2011 equ32 produced poor results
  //Fe=rhoe*E-1/2*E^2*grad(epsilon)
  //grad(epsilon)i=sum(mj/rhoj*(epsilonj-epsiloni)*grad(W))
  for(unsigned int i=0;i!=Region._PtPairList.size();++i)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         // ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         //  ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair
         )
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          if(PtiPtr!=PtjPtr)
            {
              epsiloni=PtiPtr->_eEpsilon;
              epsilonj=PtjPtr->_eEpsilon;
              
              normx=(epsilonj-epsiloni)*KnlPtr->_Wx;
              normy=(epsilonj-epsiloni)*KnlPtr->_Wy;

              // //use epsilon not smoothed
              // normx=(Region._PartList[PtjPtr->_PID-1]._eEpsilon-Region._PartList[PtiPtr->_PID-1]._eEpsilon)*KnlPtr->_Wx;
              // normy=(Region._PartList[PtjPtr->_PID-1]._eEpsilon-Region._PartList[PtiPtr->_PID-1]._eEpsilon)*KnlPtr->_Wy;
              
              PtiPtr->_Fex+=PtjPtr->_mrho*normx;
              PtiPtr->_Fey+=PtjPtr->_mrho*normy;

              if(PtjPtr!=PtiPtr&&PtjPtr->_Type==enSPHPt)
                {
                  PtjPtr->_Fex+=PtiPtr->_mrho*normx;
                  PtjPtr->_Fey+=PtiPtr->_mrho*normy;
                }
            }
        }
    }

  double ModE;
  for(unsigned int i=0;i!=Region._CalList.size() ;++i)
    {
      BasePtPtr=&Region._PtList[ID=Region._CalList[i]];

      Ex=BasePtPtr->_eEx;
      Ey=BasePtPtr->_eEy;

      ModE=Ex*Ex+Ey*Ey;
      
      BasePtPtr->_Fex*=-0.5*ModE;
      BasePtPtr->_Fey*=-0.5*ModE;

      BasePtPtr->_Fex+=BasePtPtr->_eRho*Ex;
      BasePtPtr->_Fey+=BasePtPtr->_eRho*Ey;
    }

  CKFile kfiletemp;
  kfiletemp.outTecplotEHDPLANNER(Region,Timesteps,"xxehd-case3");

  cout<<"----------------------------------------------"<<endl;
}

void CSPHEHD::Solve3(CRegion & Region, unsigned int Timesteps)
{
  CSPHPt * PtiPtr,*PtjPtr;
  CBasePt * BasePtPtr;
  CKnl * KnlPtr;
  double epsiloni, epsilonj;//epsilon of part
  double kappai, kappaj;//kappa of part
  double norm;
  double normx,normy;

  double normx1,normy1;
  double Ex,Ey;
  double Exi,Exj,Eyi,Eyj;
  unsigned int IDi,IDj;  
  unsigned int ID,ID2;

  //1.update eepsilon and ekappa
  //Lopez 2011 use an interpolation scheme to calculate transition value between the two fluid with VOF colour
  //here we try an direct interpolation using SPH, just like the interpolation of colour
 
  for (unsigned int i=0;i!=Region._PtPairList.size();i++)
    {
      // if(Region._PtPairList[i]._Type!=enSphbndptpair&&Region._PtPairList[i]._Type!=enSPHDumPtPair)
      {
        PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
        PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
        KnlPtr=&Region._KnlList[i];

        epsiloni=Region._PartList[PtiPtr->_PID-1]._eEpsilon;
        epsilonj=Region._PartList[PtjPtr->_PID-1]._eEpsilon;

        kappai=Region._PartList[PtiPtr->_PID-1]._eKappa;
        kappaj=Region._PartList[PtjPtr->_PID-1]._eKappa;

        if(PtiPtr->_Type==enSPHPt)
          {
            PtiPtr->_eEpsilon+=PtjPtr->_mrho*KnlPtr->_W*epsilonj;

            PtiPtr->_eKappa+=PtjPtr->_mrho*KnlPtr->_W*kappaj;
          }
        
        if(PtjPtr->_Type==enSPHPt)         
          {
            if(PtiPtr!=PtjPtr)
              {
                PtjPtr->_eEpsilon+=PtiPtr->_mrho*KnlPtr->_W*epsiloni;

                PtjPtr->_eKappa+=PtiPtr->_mrho*KnlPtr->_W*kappai;
              }
          }

        //for ehd dummy or ehd boundary particles
        else
          {
            PtjPtr->_eEpsilon=epsilonj;
            PtjPtr->_eKappa=kappaj;
          }
      }
    }


  //temp, for test, do not interpolate kappa and epsilon
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //temp, specify epsilon and kappa from the part, without transient region between two fluid
  // for(unsigned int i=0;i!=Region._PtList.size();++i)
  //   {
  //     BasePtPtr=&Region._PtList[i];
  //     BasePtPtr->_eEpsilon=Region._PartList[BasePtPtr->_PID-1]._eEpsilon;
  //     BasePtPtr->_eKappa=Region._PartList[BasePtPtr->_PID-1]._eKappa;
  //   }
  
  //interpolate Ex,Ey, Phie of ehd boundary particle from fluid particles
  InterpEHDBndProperty(Region, Timesteps);

  //intepolate phi of ehd dummy particles from fluid and ehd boundary particles
  InterpEHDDumPropterty(Region, Timesteps);

  //2.calulate phi n+1, construct A(phi n+1)=b, nabla((K*deltaT+epsilon)*nabla(phi n+1))=-rhoe n
  //---------------------------------------------------------------------------------------------------------------------
  //use laspack to solve linear equ.
  unsigned int PtNum=Region._CalList.size();//the number used to initializae the size of linear equation
  QMatrix LasA;//laspack amtrix A
  Vector Lasx,Lasb;//laspack vector, not std vector
  size_t Dim=PtNum;

  char namea[]="A";
  char namex[]="x";
  char nameb[]="b";
  Q_Constr(&LasA, namea, Dim , False , Rowws, Normal , True );
  V_Constr(&Lasx, namex, Dim, Normal, True);
  V_Constr(&Lasb, nameb, Dim, Normal, True);

  //set the length of LasA, after q_setlen, element will be set to 0
  //for i, some neighbor will not involved in calculation (if neighbor j is dumpt or not in calculate range)
  //so the length of row may larger than the element number, but it does not influence the result of Ax=b
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];

      Q_SetLen(&LasA, BasePtPtr->_ID2, BasePtPtr->_NumNegbor);

      Q_SetEntry(&LasA, BasePtPtr->_ID2, 0, BasePtPtr->_ID2, 0.0);//set the first one of row as i-i element
    }
  //set all components of lasvector 0
  V_SetAllCmp(&Lasb, 0.0);

  //initilize laspack matrix A and vector b, done
  //---------------------------------------------------------------------------------------------------------------------

  //3.1.2 contruct and solove Ax=b
  // construct A(phi n+1)=b, nabla((K*deltaT+epsilon)*nabla(phi n+1))=-rhoe n
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
      V_AddCmp(&Lasb, BasePtPtr->_ID2, -BasePtPtr->_eRho);
    }

  vector<unsigned int> icount;//element number of each row
  icount.resize(PtNum,1);//there is one element now, i-i, the first one of each row
  double DeltaT;
  DeltaT=Region._ControlSPH._DeltaT;

  for(unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          //as (phi i-phi j ) is calculated,  so ptiptr=ptjptr makes no sense
          if(PtiPtr!=PtjPtr)
            {
              IDi=PtiPtr->_ID2;
              IDj=PtjPtr->_ID2;

              // epsiloni=Region._PartList[PtiPtr->_PID-1]._eEpsilon;//epsilon not smoothed
              // epsilonj=Region._PartList[PtjPtr->_PID-1]._eEpsilon;

              // kappai=Region._PartList[PtiPtr->_PID-1]._eKappa;//kappa not smoothed
              // kappaj=Region._PartList[PtjPtr->_PID-1]._eKappa;

              epsiloni=PtiPtr->_eEpsilon;
              epsilonj=PtjPtr->_eEpsilon;

              kappai=PtiPtr->_eKappa;
              kappaj=PtjPtr->_eKappa;

              // norm=(kappai*DeltaT+epsiloni+kappaj*DeltaT+epsilonj)*KnlPtr->_Ww;
              norm=4*(kappai*DeltaT+epsiloni)*(kappaj*DeltaT+epsilonj)/(kappai*DeltaT+epsiloni+kappaj*DeltaT+epsilonj)*KnlPtr->_Ww;
              //norm=2*(kappai*DeltaT+epsiloni)*KnlPtr->_Ww;//test for jump between two parts
              
              Q_AddVal(&LasA, IDi, 0, PtjPtr->_mrho*norm);//the first one in row is i i

              if(PtjPtr->_Type==enSPHPt)
                {
                  icount[IDi-1]++;

                  Q_SetEntry(&LasA, IDi, icount[IDi-1]-1, IDj, -PtjPtr->_mrho*norm);

                  icount[IDj-1]++;

                  //norm=2*(kappaj*DeltaT+epsilonj)*KnlPtr->_Ww;//test for jump between two parts
                  
                  Q_AddVal(&LasA, IDj, 0, PtiPtr->_mrho*norm);//the first one in row j-j
                  Q_SetEntry(&LasA, IDj, icount[IDj-1]-1, IDi, -PtiPtr->_mrho*norm);

                }

              //ptjptr is a ehd dummy or ehd boundary particle
              else
                {
                  V_AddCmp(&Lasb, IDi, PtjPtr->_mrho*norm*PtjPtr->_ePhi);
                }
            }
        }
    }
  icount.clear();

  SetRTCAccuracy(1e-8);
  //set the initial value using the previous timestep value
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
      V_SetCmp(&Lasx, BasePtPtr->_ID2, BasePtPtr->_ePhi);
    }
  // V_SetAllCmp(&Lasx, 0.0);//initial value 0

  BiCGSTABIter(&LasA, &Lasx, &Lasb, 100, SSORPrecond,1.2);//the 4th parameter 100, is the maxiter steps
  // outputV( Lasx,  "xxbicgs-20");
  //output2(LasA, Lasx, Lasb, "xxbicgs-ax=b");

  //set φ of particle from b
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];

      BasePtPtr->_ePhi=V_GetCmp(&Lasx, BasePtPtr->_ID2);
    }

  Q_Destr(&LasA);
  V_Destr(&Lasx);
  V_Destr(&Lasb);

  //drhoe n+1=nabla*((k*deltaT)*nabla(phi n+1))
  for(unsigned int i=0;i!=Region._PtList.size();++i)
    {
      //Region._PtList[i]._eRho=0.0;
      Region._PtList[i]._deRho=0.0;
    }
  DeltaT=Region._ControlSPH._DeltaT;
  //#pragma omp parallel for
  for(unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          //drhoei=-sum(mj/rhoj*(ki+kj)*DeltaT*(phii-phij)*Ww)
          if(PtiPtr!=PtjPtr)
            {
              kappai=PtiPtr->_eKappa;
              kappaj=PtjPtr->_eKappa;

              // kappai=Region._PartList[PtiPtr->_PID-1]._eKappa;
              // kappaj=Region._PartList[PtjPtr->_PID-1]._eKappa;
              
              //norm=(kappai+kappaj)*DeltaT*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Ww;              
              //norm=2*kappai*DeltaT*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Ww;
              if(ISZERO(kappai)||ISZERO(kappaj))
                norm=0;
              else
                norm=4*kappai*kappaj*DeltaT/(kappai+kappaj)*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Ww;
              
              PtiPtr->_deRho+=PtjPtr->_mrho*norm;

              if(PtjPtr->_Type==enSPHPt)
                {
                  // norm=2*kappaj*DeltaT*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Ww;
                  
                  PtjPtr->_deRho-=PtiPtr->_mrho*norm;
                }            
            }
        }
    }

  //caculate Rhoe
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];

      // if(ISZERO(BasePtPtr->_eKappa))
      //   {
      //     BasePtPtr->_deRho=0.0;
      //   }

      BasePtPtr->_eRho+=BasePtPtr->_deRho;
    }
  
 
  //calculate E n+1=-grad(phi n+1)
  //---------------------------------------------------------------------------------------------------------------------
  //3.2 calculate E using CSPM, as the dirct scheme produce unacceptable errors
  //3.2.1 cspm coefficient matrix
  if(Region._ControlSPH._CSPMIflag2==0)
    {
      _CSPMEHD.GetCSPMGradCorctCoef2(Region);//dummy particle involved
    }

  // //3.2.2 cspm source term for E=grad(phi)
  vector<double> bx;//CSPM修正的源项
  vector<double> by;
  bx.resize(Region._CalList.size(),0.0);
  by.resize(Region._CalList.size(),0.0);
  for(unsigned int i=0;i<Region._PtPairList.size();i++)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          IDi=PtiPtr->_ID2;
          IDj=PtjPtr->_ID2;

          if(PtiPtr!=PtjPtr)
            {
              bx[IDi-1]+=PtjPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wx;
              by[IDi-1]+=PtjPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wy;

              //particle j
              if(PtjPtr->_Type==enSPHPt)
                {
                  bx[IDj-1]+=PtiPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wx;
                  by[IDj-1]+=PtiPtr->_mrho*(PtiPtr->_ePhi-PtjPtr->_ePhi)*KnlPtr->_Wy;
                }
            }
        }
    }

  //calculate CSPM corrected E=grad(phi)
  for(unsigned int i=0;i!=Region._CalList.size();++i)
    {
      BasePtPtr=&Region._PtList[Region._CalList[i]];
      ID2=BasePtPtr->_ID2;

      BasePtPtr->_eEx=-(BasePtPtr->_CSPMAxx*bx[ID2-1]+BasePtPtr->_CSPMAyx*by[ID2-1]);
      BasePtPtr->_eEy=-(BasePtPtr->_CSPMAxy*bx[ID2-1]+BasePtPtr->_CSPMAyy*by[ID2-1]);
    }

  bx.clear();
  by.clear();
  //---------------------------------------------------------------------------

  //3.3 calculate the accelerate caused by ehd force
  //use equ 11, since the derivative of E is unaccruate and formulation of Lopez 2011 equ32 produced poor results
  //Fe=rhoe*E-1/2*E^2*grad(epsilon)
  //grad(epsilon)i=sum(mj/rhoj*(epsilonj-epsiloni)*grad(W))
  for(unsigned int i=0;i!=Region._PtPairList.size();++i)
    {
      if(Region._PtPairList[i]._Type==enSPHPtPair
         // ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
         //  ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair
         )
        {
          PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
          PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
          KnlPtr=&Region._KnlList[i];

          if(PtiPtr!=PtjPtr)
            {
              epsiloni=PtiPtr->_eEpsilon;
              epsilonj=PtjPtr->_eEpsilon;
              
              normx=(epsilonj-epsiloni)*KnlPtr->_Wx;
              normy=(epsilonj-epsiloni)*KnlPtr->_Wy;

              // //use epsilon not smoothed
              // normx=(Region._PartList[PtjPtr->_PID-1]._eEpsilon-Region._PartList[PtiPtr->_PID-1]._eEpsilon)*KnlPtr->_Wx;
              // normy=(Region._PartList[PtjPtr->_PID-1]._eEpsilon-Region._PartList[PtiPtr->_PID-1]._eEpsilon)*KnlPtr->_Wy;
              
              PtiPtr->_Fex+=PtjPtr->_mrho*normx;
              PtiPtr->_Fey+=PtjPtr->_mrho*normy;

              if(PtjPtr!=PtiPtr&&PtjPtr->_Type==enSPHPt)
                {
                  PtjPtr->_Fex+=PtiPtr->_mrho*normx;
                  PtjPtr->_Fey+=PtiPtr->_mrho*normy;
                }
            }
        }
    }

  double ModE;
  for(unsigned int i=0;i!=Region._CalList.size() ;++i)
    {
      BasePtPtr=&Region._PtList[ID=Region._CalList[i]];

      Ex=BasePtPtr->_eEx;
      Ey=BasePtPtr->_eEy;

      ModE=Ex*Ex+Ey*Ey;
      
      BasePtPtr->_Fex*=-0.5*ModE;
      BasePtPtr->_Fey*=-0.5*ModE;

      BasePtPtr->_Fex+=BasePtPtr->_eRho*Ex;
      BasePtPtr->_Fey+=BasePtPtr->_eRho*Ey;
    }

  CKFile kfiletemp;
  // kfiletemp.outTecplotEHDPLANNER(Region,Timesteps,"xxehd-case3");
  if(Timesteps==0||((Timesteps+1)%Region._ControlSPH._OutputSteps)==0)
    {
       kfiletemp.outTecplotEHDDrop(Region, Timesteps+1,"xx"+Region._ControlSPH._InfileName);
       // kfiletemp.outTecplotIsoCondCylinder(Region, Timesteps+1,"xx"+Region._ControlSPH._InfileName);
      // kfiletemp.outTecplotEHDBulkRelax(Region, Timesteps+1,"xx"+Region._ControlSPH._InfileName);
    }

  cout<<"----------------------------------------------"<<endl;
}


//interpolate phi of ehd dummy particle from ehd dummy and fluid particles
void CSPHEHD::InterpEHDDumPropterty(CRegion & Region, unsigned int TimeSteps)
{
  CSPHPt * PtiPtr,*PtjPtr;
  CBasePt * PtPtr;
  CKnl * KnlPtr;
  unsigned int ID;

  //if timesteps==0, set phi of ehd dummy particles same as the ehd boundary particles
  //which has been done in model.cpp
  if(TimeSteps!=0)
    {
  
      for (unsigned int i=0;i!=Region._PtList.size();++i)
        {
          PtPtr=&Region._PtList[i];

          if(PtPtr->_Type==enEHDDumPt)
            {
              PtPtr->_CSPMCoef=0.0;

              PtPtr->_ePhi=0.0;
            }
        }

      //use Stranex's SPH to interpolate phi of ehd dummy particles
      //[sum(mj/rhoj*Wij)  sum(mj/rhoj*xij*Wij)  sum(mj/rhoj*yij*Wij)] * [f ] = [sum(mj/rhoj*fj*Wij)]
      //[sum(mj/rhoj*Wx )  sum(mj/rhoj*xij*Wx )  sum(mj/rhoj*yij*Wx )] * [fx] = [sum(mj/rhoj*fj*Wx )]
      //[sum(mj/rhoj*Wy )  sum(mj/rhoj*xij*Wy )  sum(mj/rhoj*yij*Wy )] * [fy] = [sum(mj/rhoj*fj*Wy )]

      unsigned int EHDDumNum=0;//used to count the number of ehd dummy particles
      vector<unsigned int> EHDDumID;//an idex of ehddummy particles
      EHDDumID.resize(Region._PtList.size(),0);
      for(unsigned int i=0;i!=Region._PtList.size();++i)
        {
          if(Region._PtList[i]._Type==enEHDDumPt)
            {
              EHDDumNum++;
              EHDDumID[i]=EHDDumNum;//set the i element as the index of ehd dummy particle
            }
        }

      vector<double> a11,a12,a13,a21,a22,a23,a31,a32,a33;
      vector<double> bPhi1,bPhi2,bPhi3;
      vector<double> bEx1,bEx2,bEx3;//used to calculate Ex
      vector<double> bEy1,bEy2,bEy3;//used to calculate Ey

      a11.resize(EHDDumNum,0.0);
      a12.resize(EHDDumNum,0.0);
      a13.resize(EHDDumNum,0.0);
      a21.resize(EHDDumNum,0.0);
      a22.resize(EHDDumNum,0.0);
      a23.resize(EHDDumNum,0.0);
      a31.resize(EHDDumNum,0.0);
      a32.resize(EHDDumNum,0.0);
      a33.resize(EHDDumNum,0.0);

      bPhi1.resize(EHDDumNum,0.0);
      bPhi2.resize(EHDDumNum,0.0);
      bPhi3.resize(EHDDumNum,0.0);
      bEx1.resize(EHDDumNum,0.0);
      bEx2.resize(EHDDumNum,0.0);
      bEx3.resize(EHDDumNum,0.0);
      bEy1.resize(EHDDumNum,0.0);
      bEy2.resize(EHDDumNum,0.0);
      bEy3.resize(EHDDumNum,0.0);      

      double xij,yij;
      double norm,normx,normy;
 
      for (unsigned int i=0;i<Region._PtPairList.size();i++)
        {
          if(Region._PtPairList[i]._Type==enSPHEHDDumPtPair||Region._PtPairList[i]._Type==enEHDBndDumPtPair)
            {
              PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
              PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
              KnlPtr=&Region._KnlList[i];

              xij=Region._PtPairList[i]._xij;
              yij=Region._PtPairList[i]._yij;

              norm=PtiPtr->_mrho*KnlPtr->_W;
              normx=PtiPtr->_mrho*KnlPtr->_Wx;
              normy=PtiPtr->_mrho*KnlPtr->_Wy;
              
              ID=PtjPtr->_ID;
              a11[EHDDumID[ID-1]-1]+=norm;
              a12[EHDDumID[ID-1]-1]+=xij*norm;
              a13[EHDDumID[ID-1]-1]+=yij*norm;
              a21[EHDDumID[ID-1]-1]+=normx;
              a22[EHDDumID[ID-1]-1]+=xij*normx;
              a23[EHDDumID[ID-1]-1]+=yij*normx;
              a31[EHDDumID[ID-1]-1]+=normy;
              a32[EHDDumID[ID-1]-1]+=xij*normy;
              a33[EHDDumID[ID-1]-1]+=yij*normy;

              bPhi1[EHDDumID[ID-1]-1]+=PtiPtr->_ePhi*norm;
              bPhi2[EHDDumID[ID-1]-1]+=PtiPtr->_ePhi*normx;
              bPhi3[EHDDumID[ID-1]-1]+=PtiPtr->_ePhi*normy;

              bEx1[EHDDumID[ID-1]-1]+=PtiPtr->_eEx*norm;
              bEx2[EHDDumID[ID-1]-1]+=PtiPtr->_eEx*normx;
              bEx3[EHDDumID[ID-1]-1]+=PtiPtr->_eEx*normy;
              
              bEy1[EHDDumID[ID-1]-1]+=PtiPtr->_eEy*norm;
              bEy2[EHDDumID[ID-1]-1]+=PtiPtr->_eEy*normx;
              bEy3[EHDDumID[ID-1]-1]+=PtiPtr->_eEy*normy;
    
            }
        }
      //use laspack to solve linear equ.
      QMatrix LasA;//laspack amtrix A
      Vector Lasx,Lasb;//laspack vector, not std vector
      size_t Dim=3;

      char namea[]="A";
      char namex[]="x";
      char nameb[]="b";
      Q_Constr(&LasA, namea, Dim , False , Rowws, Normal , True );
      V_Constr(&Lasx, namex, Dim, Normal, True);
      V_Constr(&Lasb, nameb, Dim, Normal, True);
      //set all components of lasvector 0
      V_SetAllCmp(&Lasb, 0.0);
      //set the length of A as Dim
      for(unsigned int i=0;i!=Dim;++i)
        {
          Q_SetLen(&LasA, i+1, Dim);
        }
      unsigned int icount=0;
      for (unsigned int i=0;i!=Region._PtList.size();++i)
          {
            PtPtr=&Region._PtList[i];

            if(PtPtr->_Type==enEHDDumPt)
              {
                Q_SetEntry(&LasA, 1, 0, 1, a11[icount]);
                Q_SetEntry(&LasA, 1, 1, 2, a12[icount]);
                Q_SetEntry(&LasA, 1, 2, 3, a13[icount]);
                Q_SetEntry(&LasA, 2, 0, 1, a21[icount]);
                Q_SetEntry(&LasA, 2, 1, 2, a22[icount]);
                Q_SetEntry(&LasA, 2, 2, 3, a23[icount]);
                Q_SetEntry(&LasA, 3, 0, 1, a31[icount]);
                Q_SetEntry(&LasA, 3, 1, 2, a32[icount]);
                Q_SetEntry(&LasA, 3, 2, 3, a33[icount]); 

                //Phi
                V_SetCmp(&Lasb, 1, bPhi1[icount]);
                V_SetCmp(&Lasb, 2, bPhi2[icount]);
                V_SetCmp(&Lasb, 3, bPhi3[icount]);
                
                SetRTCAccuracy(1e-8);
                V_SetAllCmp(&Lasx, 0.0);
 
                BiCGIter(&LasA, &Lasx, &Lasb, 100, SSORPrecond,1.2);//the 4th parameter 100, is the maxiter steps
                PtPtr->_ePhi=V_GetCmp(&Lasx, 1);
                // PtPtr->_eEx=V_GetCmp(&Lasx, 2);
                // PtPtr->_eEy=V_GetCmp(&Lasx, 3);

                //Ex
                V_SetCmp(&Lasb, 1, bEx1[icount]);
                V_SetCmp(&Lasb, 2, bEx2[icount]);
                V_SetCmp(&Lasb, 3, bEx3[icount]);

                V_SetAllCmp(&Lasx, 0.0);
 
                BiCGIter(&LasA, &Lasx, &Lasb, 100, SSORPrecond,1.2);//the 4th parameter 100, is the maxiter steps
                PtPtr->_eEx=V_GetCmp(&Lasx, 1);

                //Ey
                V_SetCmp(&Lasb, 1, bEy1[icount]);
                V_SetCmp(&Lasb, 2, bEy2[icount]);
                V_SetCmp(&Lasb, 3, bEy3[icount]);

                V_SetAllCmp(&Lasx, 0.0);
 
                BiCGIter(&LasA, &Lasx, &Lasb, 100, SSORPrecond,1.2);//the 4th parameter 100, is the maxiter steps
                PtPtr->_eEy=V_GetCmp(&Lasx, 1);

                icount++;
              }
          }
      Q_Destr(&LasA);
      V_Destr(&Lasx);
      V_Destr(&Lasb);

      a11.clear();      a12.clear();      a13.clear();
      a21.clear();      a22.clear();      a23.clear();
      a31.clear();      a32.clear();      a33.clear();
      bPhi1.clear();      bPhi2.clear();      bPhi3.clear();

      bEx1.clear();      bEx2.clear();      bEx3.clear();
      bEy1.clear();      bEy2.clear();      bEy3.clear();
    }
}

//interpolate E of ehd boundary particles from fluid particles
void CSPHEHD::InterpEHDBndProperty(CRegion & Region,unsigned int TimeSteps)
{
  CSPHPt * PtiPtr,*PtjPtr;
  CBasePt * PtPtr;
  CKnl * KnlPtr;
  unsigned int ID;

  //if timesteps==0, set phi of ehd dummy particles same as the ehd boundary particles
  //which has been done in model.cpp
  if(TimeSteps!=0)
    {  
      for (unsigned int i=0;i!=Region._PtList.size();++i)
        {
          PtPtr=&Region._PtList[i];

          if(PtPtr->_Type==enEHDBndPt)
            {
              PtPtr->_CSPMCoef=0.0;

              PtPtr->_eRho=0.0;
              
              PtPtr->_eEx=0.0;
              PtPtr->_eEy=0.0;
            }
        }

      double norm;
      for (unsigned int i=0;i<Region._PtPairList.size();i++)
        {
          if(Region._PtPairList[i]._Type==enSPHEHDBndPtPair)
            {
              PtiPtr=(CSPHPt*)Region._PtPairList[i]._PtiPtr;
              PtjPtr=(CSPHPt*)Region._PtPairList[i]._PtjPtr;
              KnlPtr=&Region._KnlList[i];

              norm=PtiPtr->_mrho*KnlPtr->_W;
              PtjPtr->_CSPMCoef+=norm;

              PtjPtr->_eEx+=norm*PtiPtr->_eEx;
              PtjPtr->_eEy+=norm*PtiPtr->_eEy;

              PtjPtr->_eRho+=norm*PtiPtr->_eRho;
            }
        }

      
      for (unsigned int i=0;i!=Region._PtList.size();++i)
        {
          PtPtr=&Region._PtList[i];

          if(PtPtr->_Type==enEHDBndPt)
            {
              if(!ISZERO(PtPtr->_CSPMCoef))
                {
                  PtPtr->_eEx/=PtPtr->_CSPMCoef;
                  PtPtr->_eEy/=PtPtr->_CSPMCoef;

                  PtPtr->_eRho/=PtPtr->_CSPMCoef;
                }
            }
        }
    }
}

void CSPHEHD::output(double ** a, double * b, int n)
{
  ofstream outfile;
  ostringstream outfilename;
  outfilename<<"Matrix-ehd"<<".dat";
  outfile.open(outfilename.str().data(),ios::out);

  if (!outfile)
    {
      cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
      return;
    }

  for (size_t i=0;i<n;i++)
    {
      for(size_t j=0;j!=n;++j)
        {
          outfile
            <<setiosflags(ios_base::scientific)
            <<setprecision(16)
            <<setw(26)<<a[j][i]<<" ";
        }

      outfile
        <<setiosflags(ios_base::scientific)
        <<setprecision(16)
        <<setw(26)<<b[i]<<" "
        <<endl;

    }


  cout<<"Output file has been written."<<endl;

  outfile.close();
}


void CSPHEHD::output2(QMatrix & a, Vector & x, Vector & b ,string outputname)
{

  ofstream outfile;
  ostringstream outfilename;
  outfilename<<outputname<<".dat";
  outfile.open(outfilename.str().data(),ios::out);

  if (!outfile)
    {
      cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
      return;
    }

    size_t n=Q_GetDim(&a);
    for (size_t i=0;i<n;i++)
      {
        for(size_t j=0;j!=n;++j)
          {
            outfile
              <<setiosflags(ios_base::scientific)
              <<setprecision(16)
              <<setw(26)<<Q_GetEl(&a, i+1, j+1)<<" ";
          }

        outfile
          <<setiosflags(ios_base::scientific)
          <<setprecision(16)
          <<setw(26)<<V_GetCmp(&x, i+1)<<" "
          <<setw(26)<<V_GetCmp(&b, i+1)<<" "
          <<endl;

      }


    cout<<"Output2 file has been written."<<endl;

    outfile.close();
  }


  void CSPHEHD::outputV(Vector & x ,string outputname)
  {

    ofstream outfile;
    ostringstream outfilename;
    outfilename<<outputname<<".dat";
    outfile.open(outfilename.str().data(),ios::out);

    if (!outfile)
      {
        cerr<<"Out put file could not be open."<<endl<<"tecplot.dat"<<endl;
        return;
      }

    size_t n=V_GetDim(&x);
    for (size_t i=0;i<n;i++)
      {
        outfile
          <<setiosflags(ios_base::scientific)
          <<setprecision(16)
          <<setw(26)<<V_GetCmp(&x, i+1)<<" "
          <<endl;
      }


    cout<<outputname<<" file has been written."<<endl;

    outfile.close();
  }
