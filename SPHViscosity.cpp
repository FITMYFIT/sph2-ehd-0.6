#include "SPHViscosity.h"

CSPHViscosity::CSPHViscosity()
{
}

CSPHViscosity::~CSPHViscosity()
{
}


void CSPHViscosity::Solve(CRegion &Region,unsigned int TimeSteps)
{
	CSPHPt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;

	double uij,vij,wij;
	double xij,yij,zij;
	double mi,mj;
	double rhoi,rhoj;
	double rhoav;
	double etai,etaj;
	double visnorm;
	unsigned int IDi,IDj;
	double wxc,wyc,wzc;//修正核函数梯度
	double MAXX,MINX;//用于周期性边界条件

	if (Region._ControlSPH._PerdBnd==1)
	{
		MAXX=Region._ControlSPH._PerdBndMaxX;
		MINX=Region._ControlSPH._PerdBndMinX;
	}


	SolveEta(Region,TimeSteps);//牛顿流体时:Eta=Mu=k;幂律流：Eta=k*gamma^(n-1)

	for (size_t i=0;i<Region._PtPairList.size();i++)
	{
		PtiPtr=(CSPHPt *)Region._PtPairList[i]._PtiPtr;
		PtjPtr=(CSPHPt *)Region._PtPairList[i]._PtjPtr;
		KnlPtr=&Region._KnlList[i];

		if(Region._PtPairList[i]._Type==enSPHPtPair)
		{
			if(PtiPtr!=PtjPtr)
			{
				IDi=PtiPtr->_ID2;
				IDj=PtjPtr->_ID2;

				mi=PtiPtr->_m;
				mj=PtjPtr->_m;

				rhoi=PtiPtr->_rho;
				rhoj=PtjPtr->_rho;

				//幂律型表观粘度
				etai=PtiPtr->_VisEta;
				etaj=PtjPtr->_VisEta;

				xij=Region._PtPairList[i]._xij;
				yij=Region._PtPairList[i]._yij;

				if(Region._ControlSPH._RunMod==1)//显式算法
          {
            uij=PtiPtr->_u-PtjPtr->_u;
            vij=PtiPtr->_v-PtjPtr->_v;

            //particle i 粘性项
            visnorm=mj*(etai+etaj)/(rhoi*rhoj)*KnlPtr->_Ww;//Morris
            //visnorm=4*mj*(etai+etaj)/(rhoi*rhoi+rhoj*rhoj)*KnlPtr->_Ww;//X J Fan
            //visnorm=4*mj/(rhoi*rhoj)*(etai*etaj)/(etai+etaj)*KnlPtr->_Ww;//P W Cleary热传导类似
            //visnorm=mj*(etai+etaj)/(rhoi*rhoj)/(xij*xij+yij*yij+zij*zij)*(xij*KnlPtr->_Wxci+yij*KnlPtr->_Wyci+zij*KnlPtr->_Wzci);//SPH-4 P134 式4

            PtiPtr->_du+=uij*visnorm;
            PtiPtr->_dv+=vij*visnorm;

            //particle j 粘性项
            visnorm=mi*(etai+etaj)/(rhoi*rhoj)*KnlPtr->_Ww;//Morris
            //visnorm=4*mi*(etai+etaj)/(rhoi*rhoi+rhoj*rhoj)*KnlPtr->_Ww;//X J Fan
            //visnorm=4*mi/(rhoi*rhoj)*(etai*etaj)/(etai+etaj)*KnlPtr->_Ww;//P W Cleary热传导类似
            //visnorm=mi*(etai+etaj)/(rhoi*rhoj)/(xij*xij+yij*yij+zij*zij)*(xij*KnlPtr->_Wxcj+yij*KnlPtr->_Wycj+zij*KnlPtr->_Wzcj);//SPH-4 P134 式4

            PtjPtr->_du-=uij*visnorm;
            PtjPtr->_dv-=vij*visnorm;
          }

				else//隐式算法
				{
					CMatrix TempM;

					//particle i 粘性项
					//----------------------------------------------------
					visnorm=mj*(etai+etaj)/(rhoi*rhoj)*KnlPtr->_Ww;//Morris
					//visnorm=4*mj*(etai+etaj)/(rhoi*rhoi+rhoj*rhoj)*KnlPtr->_Ww;//X J Fan
					//visnorm=4*mj/(rhoi*rhoj)*(etai*etaj)/(etai+etaj)*KnlPtr->_Ww;//P W Cleary热传导类似

					////两粒子光滑长度不一致时的计算式R Vacondio Equ.7(略有改动)
					//rhoav=0.5*(rhoi+rhoj);
					//visnorm=mj*(etai+etaj)/(rhoav*rhoav)*0.5*(KnlPtr->_Ww+KnlPtr->_Wwj);


					//_GetKnlListV.GetCrctKnlGrad(Region,PtiPtr,KnlPtr,wxc,wyc,wzc);//求出修正核函数
					//visnorm=mj*(etai+etaj)/(rhoi*rhoj)/(xij*xij+yij*yij+zij*zij)*(xij*wxc+yij*wyc+zij*wzc);//SPH-4 P134 式4

					Region._SPHAm[IDi-1]._Ele-=visnorm;

					TempM._Ele=visnorm;
					TempM._RowID=IDi;
					TempM._ColID=IDj;
					Region._SPHAm.push_back(TempM);
					//----------------------------------------------------

					//particle j 粘性项
					//----------------------------------------------------
					visnorm=mi*(etai+etaj)/(rhoi*rhoj)*KnlPtr->_Ww;//Morris
					//visnorm=4*mi*(etai+etaj)/(rhoi*rhoi+rhoj*rhoj)*KnlPtr->_Ww;//X J Fan
					//visnorm=4*mi/(rhoi*rhoj)*(etai*etaj)/(etai+etaj)*KnlPtr->_Ww;//P W Cleary热传导类似

					////两粒子光滑长度不一致时的计算式R Vacondio Equ.7(略有改动)
					//rhoav=0.5*(rhoi+rhoj);
					//visnorm=mi*(etai+etaj)/(rhoav*rhoav)*0.5*(KnlPtr->_Ww+KnlPtr->_Wwj);

					//_GetKnlListV.GetCrctKnlGrad(Region,PtjPtr,KnlPtr,wxc,wyc,wzc);//求出修正核函数
					//visnorm=mi*(etai+etaj)/(rhoi*rhoj)/(xij*xij+yij*yij+zij*zij)*(xij*wxc+yij*wyc+zij*wzc);//SPH-4 P134 式4

					Region._SPHAm[IDj-1]._Ele-=visnorm;

					TempM._Ele=visnorm;
					TempM._RowID=IDj;
					TempM._ColID=IDi;
					Region._SPHAm.push_back(TempM);
					//----------------------------------------------------
				}
			}
		}


		if(Region._PtPairList[i]._Type==enSPHNULLPtPair
       ||Region._PtPairList[i]._Type==enSPHBndPtPair
       ||Region._PtPairList[i]._Type==enSPHDumPtPair
       ||Region._PtPairList[i]._Type==enSPHEHDBndPtPair
       ||Region._PtPairList[i]._Type==enSPHEHDDumPtPair)
		{
			IDi=PtiPtr->_ID2;
			IDj=PtjPtr->_ID2;

			mi=PtiPtr->_m;
			mj=PtjPtr->_m;

			rhoi=PtiPtr->_rho;
			rhoj=PtjPtr->_rho;

			//幂律型表观粘度
			etai=PtiPtr->_VisEta;
			etaj=PtjPtr->_VisEta;

      xij=Region._PtPairList[i]._xij;
      yij=Region._PtPairList[i]._yij;
      
      if(Region._ControlSPH._RunMod==1)//显式算法
        {
          uij=PtiPtr->_u-PtjPtr->_u;
          vij=PtiPtr->_v-PtjPtr->_v;

          //particle i 粘性项
          visnorm=mj*(etai+etaj)/(rhoi*rhoj)*KnlPtr->_Ww;//Morris
          ////visnorm=4*mj*(etai+etaj)/(rhoi*rhoi+rhoj*rhoj)*KnlPtr->_Ww;//X J Fan
          ////visnorm=4*mj/(rhoi*rhoj)*(etai*etaj)/(etai+etaj)*KnlPtr->_Ww;//P W Cleary热传导类似
          ////visnorm=mj*(etai+etaj)/(rhoi*rhoj)/(xij*xij+yij*yij+zij*zij)*(xij*KnlPtr->_Wxci+yij*KnlPtr->_Wyci+zij*KnlPtr->_Wzci);//SPH-4 P134 式4

          PtiPtr->_du+=uij*visnorm;
          PtiPtr->_dv+=vij*visnorm;
        }

      else//隐式算法
        {
          CMatrix TempM;

          //particle i 粘性项
          //----------------------------------------------------
          visnorm=mj*(etai+etaj)/(rhoi*rhoj)*KnlPtr->_Ww;//Morris
          //visnorm=4*mj*(etai+etaj)/(rhoi*rhoi+rhoj*rhoj)*KnlPtr->_Ww;//X J Fan
          //visnorm=4*mj/(rhoi*rhoj)*(etai*etaj)/(etai+etaj)*KnlPtr->_Ww;//P W Cleary热传导类似

          ////两粒子光滑长度不一致时的计算式R Vacondio Equ.7(略有改动)
          //rhoav=0.5*(rhoi+rhoj);
          //visnorm=mj*(etai+etaj)/(rhoav*rhoav)*0.5*(KnlPtr->_Ww+KnlPtr->_Wwj);


          //_GetKnlListV.GetCrctKnlGrad(Region,PtiPtr,KnlPtr,wxc,wyc,wzc);//求出修正核函数
          //visnorm=mj*(etai+etaj)/(rhoi*rhoj)/(xij*xij+yij*yij+zij*zij)*(xij*wxc+yij*wyc+zij*wzc);//SPH-4 P134 式4

          Region._SPHAm[IDi-1]._Ele-=visnorm;

          Region._SPHbu[IDi-1]-=PtjPtr->_u*visnorm;
          Region._SPHbu[IDi-1]-=PtjPtr->_v*visnorm;
          //----------------------------------------------------
        }
		}
	}
}

void CSPHViscosity::SolveEta(CRegion &Region,unsigned int TimeSteps)
{
	double Ww;
	double xji,yji,zji;
	double uji,vji,wji;
	double mi,rhoi;
	double mj,rhoj;
	double normx,normy,normz;
	double wxc,wyc,wzc;//修正核函数梯度

	unsigned int IDi,IDj;
	unsigned int ID,PID;	

	CSPHPt * PtiPtr,* PtjPtr;
	CKnl * KnlPtr;

	if(Region._ControlSPH._SPHPV==0&&TimeSteps==Region._ControlSPH._StartStep)//_SPHPV==0代表需要计算牛顿型粘性,只在计算开始的第一个时间步赋值
	{
		for(unsigned int i=0;i!=Region._PtList.size();++i)
		{
			PID=Region._PtList[i]._PID;

			Region._PtList[i]._VisEta=Region._PartList[PID-1]._VisK;
		}
	}


	if(Region._ControlSPH._SPHPV==1)//_SPHPV==1代表需要计算幂律型粘性
	{
		//每个速度分量在三个方向的梯度，计算gamma的时候用
		_Gradux.resize(Region._StatDataList._InvolvedFluidNum,0.0);
		_Graduy.resize(Region._StatDataList._InvolvedFluidNum,0.0);
		_Gradvx.resize(Region._StatDataList._InvolvedFluidNum,0.0);
		_Gradvy.resize(Region._StatDataList._InvolvedFluidNum,0.0);


		//以下部分用的是sum(mj/rhoj*uji*Wx)
		//1.1 计算流体粒子的速度梯度，grad(u),grad(v)
		for(unsigned int i=0;i<Region._PtPairList.size();i++)
		{
			PtiPtr=(CSPHPt *)Region._PtPairList[i]._PtiPtr;
			PtjPtr=(CSPHPt *)Region._PtPairList[i]._PtjPtr;
			KnlPtr=&Region._KnlList[i];
			
			IDi=PtiPtr->_ID2;
			IDj=PtjPtr->_ID2;

			rhoi=PtiPtr->_rho;
			rhoj=PtjPtr->_rho;

			mi=PtiPtr->_m;
			mj=PtjPtr->_m;

			uji=PtjPtr->_u-PtiPtr->_u;
			vji=PtjPtr->_v-PtiPtr->_v;
			
			xji=-Region._PtPairList[i]._xij;
			yji=-Region._PtPairList[i]._yij;

			if(Region._PtPairList[i]._Type==enSPHPtPair)
			{	
				//--------------------------------
				//equ. graduy=sum(mj/rhoj*uji*Wy)
				normx=mj/rhoj*KnlPtr->_Wx;
				normy=mj/rhoj*KnlPtr->_Wy;

				//_GetKnlListV.GetCrctKnlGrad(Region,PtiPtr,KnlPtr,wxc,wyc,wzc);//求出修正核函数

				//normx=mj/rhoj*wxc;
				//normy=mj/rhoj*wyc;
				//normz=mj/rhoj*wzc;

				//particle i
				_Gradux[IDi-1]+=normx*uji;
				_Graduy[IDi-1]+=normy*uji;
				_Gradvx[IDi-1]+=normx*vji;
				_Gradvy[IDi-1]+=normy*vji;
							
				//particle j
				normx=mi/rhoi*KnlPtr->_Wx;
				normy=mi/rhoi*KnlPtr->_Wy;
			
				//_GetKnlListV.GetCrctKnlGrad(Region,PtjPtr,KnlPtr,wxc,wyc,wzc);//求出修正核函数

				//normx=mi/rhoi*wxc;
				//normy=mi/rhoi*wyc;
				//normz=mi/rhoi*wzc;

				_Gradux[IDj-1]+=normx*uji;
				_Graduy[IDj-1]+=normy*uji;
				_Gradvx[IDj-1]+=normx*vji;
				_Gradvy[IDj-1]+=normy*vji;
				//----------------------------------
			}

			if(Region._PtPairList[i]._Type==enSPHNULLPtPair)
			{		
				//--------------------------------
				//equ. graduy=sum(mj/rhoj*uji*Wy)
				normx=mj/rhoj*KnlPtr->_Wx;
				normy=mj/rhoj*KnlPtr->_Wy;

				//_GetKnlListV.GetCrctKnlGrad(Region,PtiPtr,KnlPtr,wxc,wyc,wzc);//求出修正核函数

				//normx=mj/rhoj*wxc;
				//normy=mj/rhoj*wyc;
				//normz=mj/rhoj*wzc;

				//particle i
				_Gradux[IDi-1]+=normx*uji;
				_Graduy[IDi-1]+=normy*uji;
				_Gradvx[IDi-1]+=normx*vji;
				_Gradvy[IDi-1]+=normy*vji;
							
				////particle j
				//normx=mi/rhoi*KnlPtr->_Wx;
				//normy=mi/rhoi*KnlPtr->_Wy;
			
				////_GetKnlListV.GetCrctKnlGrad(Region,PtjPtr,KnlPtr,wxc,wyc,wzc);//求出修正核函数

				////normx=mi/rhoi*wxc;
				////normy=mi/rhoi*wyc;
				////normz=mi/rhoi*wzc;

				//_Gradux[IDj-1]+=normx*uji;
				//_Graduy[IDj-1]+=normy*uji;
				//_Gradvx[IDj-1]+=normx*vji;
				//_Gradvy[IDj-1]+=normy*vji;
				//----------------------------------
			}

			//if(region._PtPairList[i]._Type==enSPHBndPtPair)
			//{			
			//}
		}


		double PowerLawk,PowerLawn;
		
		//1.2 计算流体粒子的剪切速率及表观粘度
		for(unsigned int i=0;i<Region._CalList.size();i++)
		{
			ID=Region._CalList[i];

			if(Region._PtList[ID]._Type==enSPHPt)
			{
				IDi=Region._PtList[ID]._ID2;

				//2.1 计算每个粒子的剪切速率
				Region._PtList[ID]._VisGamma=sqrt(2*_Gradux[IDi-1]*_Gradux[IDi-1]+2*_Gradvy[IDi-1]*_Gradvy[IDi-1]
																		  	+(_Gradvx[IDi-1]+_Graduy[IDi-1])*(_Gradvx[IDi-1]+_Graduy[IDi-1]));
			
				//2.2 计算每个粒子的表观粘度

				PID=Region._PtList[ID]._PID;

				PowerLawk=Region._PartList[PID-1]._VisK;
				PowerLawn=Region._PartList[PID-1]._VisN;
																		 
				if(Region._PtList[ID]._VisGamma==0)
					Region._PtList[ID]._VisEta=10000;
				else
					Region._PtList[ID]._VisEta=PowerLawk*pow(Region._PtList[ID]._VisGamma,PowerLawn-1);

				if(Region._PtList[ID]._VisEta>10000)
					Region._PtList[ID]._VisEta=10000;
			}
		}			

		////直接重构表观粘度
		//if(Region._ControlSPH._CSPMIflag1==0)
		//{
		//	_CSPMV.GetCSPMFunCorctCoef(Region);
		//}

		////插值得到表观粘度
		//std::vector<double> TempEta;
		//TempEta.resize(Region._CalList.size(),0.0);

		//for(unsigned int i=0;i!=Region._PtPairList.size();++i)
		//{
		//	PtiPtr=(CSPHPt *)Region._PtPairList[i]._PtiPtr;
		//	PtjPtr=(CSPHPt *)Region._PtPairList[i]._PtjPtr;
		//	KnlPtr=&Region._KnlList[i];

		//	IDi=PtiPtr->_ID2;
		//	IDj=PtjPtr->_ID2;

		//	TempEta[IDi-1]+=PtjPtr->_m/PtjPtr->_rho*KnlPtr->_W*PtiPtr->_VisEta;

		//	if(PtiPtr!=PtjPtr)
		//	{
		//		TempEta[IDj-1]+=PtiPtr->_m/PtiPtr->_rho*KnlPtr->_W*PtjPtr->_VisEta;
		//	}
		//}

		//for(unsigned int i=0;i!=Region._PtList.size();++i)
		//{
		//	if(Region._PtList[i]._Type==enSPHPt&&Region._PtList[i]._Iflag==1)
		//	{
		//		ID=Region._PtList[i]._ID2;

		//		Region._PtList[i]._VisEta=TempEta[ID-1]/Region._PtList[i]._CSPMCoef;
		//	}
		//}

		_Gradux.clear();
		_Graduy.clear();
		_Gradvx.clear();
		_Gradvy.clear();
	}
}
