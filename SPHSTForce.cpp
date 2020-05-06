#include "SPHSTForce.h"

CSPHSTForce::CSPHSTForce()
{
}

CSPHSTForce::~CSPHSTForce()
{
}

void CSPHSTForce::Solve(CRegion &Region)
{
	CBasePt * PtiPtr,* PtjPtr;
	CBasePt * BasePtr;
	CKnl * KnlPtr;

	double mi,mj;
	double rhoi,rhoj;
	double Nnxx,Nnyy;
	unsigned int C0i,C0j;
	unsigned int IDi,IDj;
	unsigned int ID,ID2;
	unsigned int PID;
	double Axx,Ayx,Axy,Ayy;
	double xij,yij;
	double normx,normy;

	vector<double> bx;//CSPM������Դ��
	vector<double> by;
	vector<unsigned int> N;//�����ж��Ƿ�Լ��㷨����Ч

	bx.resize(Region._CalList.size(),0.0);
	by.resize(Region._CalList.size(),0.0);
	N.resize(Region._CalList.size(),0);

	//1.��ֵ�õ�ɫֵ
	for(unsigned int i=0;i!=Region._PtPairList.size();++i)
	{
		PtiPtr=Region._PtPairList[i]._PtiPtr;
		PtjPtr=Region._PtPairList[i]._PtjPtr;
		KnlPtr=&Region._KnlList[i];

		C0i=Region._PartList[PtiPtr->_PID-1]._C0;
		C0j=Region._PartList[PtjPtr->_PID-1]._C0;

		PtiPtr->_C+=PtjPtr->_m/PtjPtr->_rho*KnlPtr->_W*C0j;

		if(PtiPtr!=PtjPtr)
		{
			PtjPtr->_C+=PtiPtr->_m/PtiPtr->_rho*KnlPtr->_W*C0i;
		}
	}

	//2 ���㷨�������ķ���
	//2.1����CSPM����ϵ��
	if(Region._ControlSPH._CSPMIflag2==0)
	{
		_CSPMST.GetCSPMGradCorctCoef(Region);
	}

	//2.2 ����CSPM���������Դ��
	for(unsigned int i=0;i!=Region._PtPairList.size();++i)
	{
		if(Region._PtPairList[i]._Type==enSPHPtPair)
		{
			PtiPtr=Region._PtPairList[i]._PtiPtr;
			PtjPtr=Region._PtPairList[i]._PtjPtr;
			KnlPtr=&Region._KnlList[i];

			IDi=PtiPtr->_ID2;
			IDj=PtjPtr->_ID2;

			mi=PtiPtr->_m;
			mj=PtjPtr->_m;

			rhoi=PtiPtr->_rho;
			rhoj=PtjPtr->_rho;

			normx=mj/rhoj*KnlPtr->_Wx;
			normy=mj/rhoj*KnlPtr->_Wy;

			bx[IDi-1]+=normx*(PtiPtr->_C-PtjPtr->_C);
			by[IDi-1]+=normy*(PtiPtr->_C-PtjPtr->_C);

			if(PtiPtr!=PtjPtr)
			{
				normx=mi/rhoi*KnlPtr->_Wx;
				normy=mi/rhoi*KnlPtr->_Wy;

				bx[IDj-1]-=normx*(PtjPtr->_C-PtiPtr->_C);
				by[IDj-1]-=normy*(PtjPtr->_C-PtiPtr->_C);
			}
		}
	}

	//2.3 ����CSPM������ķ������򻯷���
	double mod;
	double xi;
	for(unsigned int i=0;i!=Region._CalList.size();++i)
	{
		ID=Region._CalList[i];

		BasePtr=&Region._PtList[ID];

		ID2=BasePtr->_ID2;//ID2ʵ���ϵ���i+1

		Axx=BasePtr->_CSPMAxx;
		Ayx=BasePtr->_CSPMAyx;
		Axy=BasePtr->_CSPMAxy;
		Ayy=BasePtr->_CSPMAyy;

		_OperatorST.Reve2ndMat(Axx,Ayx,Axy,Ayy,&Axx,&Ayx,&Axy,&Ayy);

		BasePtr->_nx=Axx*bx[ID2-1]+Ayx*by[ID2-1];
		BasePtr->_ny=Axy*bx[ID2-1]+Ayy*by[ID2-1];

		//���Renormalization Matrix B����С����ֵ S Marrone2010Equ.��1��
		Axx=-Axx; Ayx=-Ayx;//B �����Ԫ����CSPMϵ������Ԫ�ص��෴��
		Axy=-Axy; Ayy=-Ayy;

		_OperatorST.Reve2ndMat(Axx,Ayx,Axy,Ayy,&Axx,&Ayx,&Axy,&Ayy);

		BasePtr->_CSPMminLambda=_OperatorST.MinEig2ndMat(Axx,Ayx,Axy,Ayy);

		//���򻯷���
		mod=sqrt(BasePtr->_nx*BasePtr->_nx+BasePtr->_ny*BasePtr->_ny);

		xi=0.01/BasePtr->_h;

		if(mod>xi)
		{
			N[ID2-1]=1;
			BasePtr->_Nnx=BasePtr->_nx/mod;
			BasePtr->_Nny=BasePtr->_ny/mod;
		}

		else
		{
			N[ID2-1]=0;
			BasePtr->_Nnx=0.0;
			BasePtr->_Nny=0.0;
		}
	}

	bx.clear();
	by.clear();

	//3 ��������
	//3.1 ����CSPM�������ʵ�Դ�ϵ������ԭCSPM�����ϵ�����������ˣ����ӶԲ�ͬ��
	vector<double> bxx;//CSPM����Nnx��Դ��
	vector<double> bxy;
	vector<double> byx;//CSPM����Nny��Դ��
	vector<double> byy;

	bxx.resize(Region._CalList.size(),0.0);
	bxy.resize(Region._CalList.size(),0.0);
	byx.resize(Region._CalList.size(),0.0);
	byy.resize(Region._CalList.size(),0.0);
	for(unsigned int i=0;i!=Region._CalList.size();++i)
	{
		ID=Region._CalList[i];

		BasePtr=&Region._PtList[ID];

		BasePtr->_CSPMAxx=0.0;
		BasePtr->_CSPMAyx=0.0;
		BasePtr->_CSPMAxy=0.0;
		BasePtr->_CSPMAyy=0.0;
	}

	for(unsigned int i=0;i!=Region._PtPairList2.size();++i)
	{
		if(Region._PtPairList2[i]._Type==enSPHPtPair)
		{
			PtiPtr=Region._PtPairList2[i]._PtiPtr;
			PtjPtr=Region._PtPairList2[i]._PtjPtr;
			KnlPtr=&Region._KnlList2[i];

			IDi=PtiPtr->_ID2;
			IDj=PtjPtr->_ID2;

			if(N[IDi-1]==1&&N[IDj-1]==1)
			{
				mi=PtiPtr->_m;
				mj=PtjPtr->_m;

				rhoi=PtiPtr->_rho;
				rhoj=PtjPtr->_rho;

        normx=mj/rhoj*KnlPtr->_Wx;
        normy=mj/rhoj*KnlPtr->_Wy;

        PtiPtr->_CSPMAxx+=(PtiPtr->_x-PtjPtr->_x)*normx;
        PtiPtr->_CSPMAyx+=(PtiPtr->_y-PtjPtr->_y)*normx;
        PtiPtr->_CSPMAxy+=(PtiPtr->_x-PtjPtr->_x)*normy;
        PtiPtr->_CSPMAyy+=(PtiPtr->_y-PtjPtr->_y)*normy;

        bxx[IDi-1]+=(PtiPtr->_Nnx-PtjPtr->_Nnx)*normx;
        bxy[IDi-1]+=(PtiPtr->_Nnx-PtjPtr->_Nnx)*normy;

        byx[IDi-1]+=(PtiPtr->_Nny-PtjPtr->_Nny)*normx;
        byy[IDi-1]+=(PtiPtr->_Nny-PtjPtr->_Nny)*normy;

        if(PtiPtr!=PtjPtr)
          {
            normx=mi/rhoi*KnlPtr->_Wx;
            normy=mi/rhoi*KnlPtr->_Wy;

            PtjPtr->_CSPMAxx-=(PtjPtr->_x-PtiPtr->_x)*normx;
            PtjPtr->_CSPMAyx-=(PtjPtr->_y-PtiPtr->_y)*normx;
            PtjPtr->_CSPMAxy-=(PtjPtr->_x-PtiPtr->_x)*normy;
            PtjPtr->_CSPMAyy-=(PtjPtr->_y-PtiPtr->_y)*normy;

            bxx[IDj-1]-=(PtjPtr->_Nnx-PtiPtr->_Nnx)*normx;
            bxy[IDj-1]-=(PtjPtr->_Nnx-PtiPtr->_Nnx)*normy;

            byx[IDj-1]-=(PtjPtr->_Nny-PtiPtr->_Nny)*normx;
            byy[IDj-1]-=(PtjPtr->_Nny-PtiPtr->_Nny)*normy;
          }
			}
		}
	}

	//3.2 ����CSPM�����������&�����������ٶ�
	for(unsigned int i=0;i!=Region._CalList.size();++i)
	{
		ID=Region._CalList[i];

		BasePtr=&Region._PtList[ID];

		ID2=BasePtr->_ID2;//ID2ʵ���ϵ���i+1

		Axx=BasePtr->_CSPMAxx;
		Ayx=BasePtr->_CSPMAyx;
		Axy=BasePtr->_CSPMAxy;
		Ayy=BasePtr->_CSPMAyy;

		_OperatorST.Reve2ndMat(Axx,Ayx,Axy,Ayy,&Axx,&Ayx,&Axy,&Ayy);

		Nnxx=Axx*bxx[ID2-1]+Ayx*bxy[ID2-1];
		Nnyy=Axy*byx[ID2-1]+Ayy*byy[ID2-1];

		BasePtr->_Curvature=-(Nnxx+Nnyy);//���ʼ������


		//�����ɱ�����ɵļ��ٶ�
		double sigma;
		PID=BasePtr->_PID;
		sigma=Region._PartList[PID-1]._STSigma;

		normx=sigma/BasePtr->_rho*BasePtr->_Curvature*BasePtr->_nx;
		normy=sigma/BasePtr->_rho*BasePtr->_Curvature*BasePtr->_ny;

		if(Region._ControlSPH._RunMod==1)//��ʽ�㷨
		{
			BasePtr->_du+=normx;
			BasePtr->_dv+=normy;
		}

		else//��ʽ�㷨
		{
			Region._SPHbu[ID2-1]+=normx;
			Region._SPHbv[ID2-1]+=normy;
		}
	}

	bxx.clear();
	bxy.clear();
	byx.clear();
	byy.clear();

	N.clear();
}
