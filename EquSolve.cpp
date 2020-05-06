#include "EquSolve.h"

CEquSolve::CEquSolve()
{
}

CEquSolve::~CEquSolve()
{
}

void CEquSolve::BICGSolve(const std::vector<CMatrix> &A, std::vector<double> &x, const std::vector<double> &b, unsigned int MaxIter, double RTC)
{
	int Iter;
  double Alpha, Beta, Rho, RhoOld = 0.0;
  double bNorm;
  size_t Dim;

	Dim = b.size();

	std::vector<double> r;
	std::vector<double> r_;
	std::vector<double> p;
	std::vector<double> p_;
	std::vector<double> q;	

	r.resize(Dim,0);
	r_.resize(Dim,0);
	p.resize(Dim,0);
	p_.resize(Dim,0);
	q.resize(Dim,0);

	std::vector<double> vres1;
	std::vector<double> vres2;
	
	vres1.resize(Dim,0);
	vres2.resize(Dim,0);

	bNorm=_Operator.SumSquare(b);//元素的平方和再开平方

	Iter=0;

  /* r = b - A * x(i) */
	if(!IsZero(_Operator.SumFabs(x)/Dim))//if (!IsZero(l1Norm_V(x) / Dim))
	{
		_Operator.Multi_MV(A,x,vres1);
		_Operator.Subs_VV(b,vres1,vres2);
		r=vres2;
	} 
  else 
	{
		r=b;
  }

  /* plain BiCG (z = r, z_ = r_) */
	r_=r;

	while((!RTCCtrl(_Operator.SumSquare(r), bNorm,RTC))&& Iter < MaxIter)
	{
		Iter++;

		Rho = _Operator.Multi_VV(r, r_);

		if (IsZero(Rho)) 
		{
		  cout<<"Error when solve matrix equation"<<endl;
		  cout<<"Error:Rho=0"<<endl;
		  break;
		}

		if (Iter == 1) 
		{
			p=r;
			p_=r_;
		}

		else 
		{
			Beta = Rho / RhoOld;

			_Operator.Multi_SV(Beta, p,vres1);
			_Operator.Addi_VV(r,vres1,vres2);
			p=vres2;
		
			_Operator.Multi_SV(Beta,p_,vres1);
			_Operator.Addi_VV(r_,vres1,vres2);
			p_=vres2;
		}

		_Operator.Multi_MV(A,p,q);

		Alpha = Rho/_Operator.Multi_VV(p_,q);

		_Operator.Multi_SV(Alpha,p,vres1);
		_Operator.AddAsign_VV(x,vres1);

		_Operator.Multi_SV(Alpha,q,vres1);
		_Operator.SubAsign_VV(r,vres1);

		_Operator.Multi_MtV(A,p_,vres1);
		_Operator.Multi_SV(Alpha,vres1,vres2);
		_Operator.SubAsign_VV(r_,vres2);

		RhoOld = Rho;
	}
	
	r.clear();
	r_.clear();
	p.clear();
	p_.clear();
	q.clear();

	vres1.clear();
	vres2.clear();
}

unsigned int CEquSolve::RTCCtrl(double rNorm,double bNorm,double RTC)//精度控制，达到精度后退出
{
	unsigned int Result;
      
	if (rNorm < RTC*bNorm || (IsZero(bNorm) && IsOne(1.0 + rNorm)))
      Result = 1;
  else
      Result = 0;


    return(Result);
}
