#include "Operator.h"

COperator::COperator()
{
}

COperator::~COperator()
{
}

double COperator::SumSquare(const std::vector<double> &v)
{
	double sum;
  double SRes;
	double cmp;

	sum=0.0;

	for(unsigned int i=0;i<v.size();i++)
	{
		cmp=v[i];
		sum+=cmp*cmp;
	}

	SRes=sqrt(sum);

	return(SRes);
}

double COperator::SumFabs(const std::vector<double> &v)//数组的元素的绝对值之和
{
	unsigned int i;
	double cmp;
  double SRes;

	SRes=0.0;
	for(i=0;i<v.size();i++)
	{
		cmp=v[i];

		SRes+=fabs(cmp);
	}	
	
	return(SRes);
}

double COperator::Multi_VV(const std::vector<double> &v1,const std::vector<double> &v2)//3.res=v1*v2
{
	unsigned int i;
	double SRes;
	double cmp1,cmp2;

	if(v1.size()==v2.size())
	{
		SRes=0.0;

		for(i=0;i<v1.size();i++)
		{
			cmp1=v1[i];
			cmp2=v2[i];

			SRes+=cmp1*cmp2;
		}

		return(SRes);
	}

	else
	{
		cout<<"Error:不同维数的向量相乘."<<endl;
		return(0);
	}
}

void COperator::AddAsign_VV(std::vector<double> &v1,const std::vector<double> &v2)//4.v1+=v2
{
	unsigned int i;

	if(v1.size()==v2.size())
	{
		for(i=0;i<v1.size();i++)
		{
			v1[i]+=v2[i];
		}
	}

	else
	{
		cout<<"Error:不同维数的向量相加."<<endl;
	}
}

void COperator::SubAsign_VV(std::vector<double> &v1,const std::vector<double> &v2)//5.v1-=v2
{
	unsigned int i;

	if(v1.size()==v2.size())
	{
		for(i=0;i<v1.size();i++)
		{
			v1[i]-=v2[i];
		}
	}

	else
	{
		cout<<"Error:不同维数的向量相减."<<endl;
	}
}

void COperator::Subs_VV(const std::vector<double> &v1,const std::vector<double> &v2,std::vector<double> &vres)//6.res=v1-v2
{
	if(v1.size()==v2.size())
	{
		if(v1.size()==vres.size())
		{
			for(unsigned int i=0;i!=v1.size();++i)
			{
				vres[i]=v1[i]-v2[i];
			}
		}

		else
		{
			cout<<"Error:结果向量的维数与相减向量的维数不同."<<endl;
		}
	}

	else
	{
		cout<<"Error:不同维数的向量相减."<<endl;
	}
}

void COperator::Addi_VV(const std::vector<double> &v1,const std::vector<double> &v2,std::vector<double> &vres)//7.vres=v1+v2
{
	if(v1.size()==v2.size())
	{
		if(v1.size()==vres.size())
		{
			for(unsigned int i=0;i!=v1.size();++i)
			{
				vres[i]=v1[i]+v2[i];
			}
		}

		else
		{
			cout<<"Error:结果向量的维数与相加向量的维数不同."<<endl;
		}
	}

	else
	{
		cout<<"Error:不同维数的向量相加."<<endl;
	}
}

void COperator::Multi_SV(double S,const std::vector<double> &v,std::vector<double> &vres)//8.vres=S*v
{
	if(v.size()==vres.size())
	{
		for(unsigned int i=0;i!=v.size();++i)
		{
			vres[i]=S*v[i];
		}
	}

	else
	{
		cout<<"Error:结果向量与计算向量维数不同."<<endl;
	}
}

void COperator::Multi_MV(const std::vector<CMatrix> &A,const std::vector<double> &b,std::vector<double> &vres)//9.SRes=A*b
{
	unsigned int i;

	unsigned int rowid,columnid;

	if(b.size()==vres.size())
	{
		for(i=0;i<vres.size();i++)
		{
			vres[i]=0.0;
		}

		for(i=0;i<A.size();i++)
		{
			rowid=A[i]._RowID;
			columnid=A[i]._ColID;

			vres[rowid-1]+=A[i]._Ele*b[columnid-1];
		}
	}
		
	else
	{
		cout<<"Error:结果向量的维数与相乘的矩阵和向量不同."<<endl;
	}
}


void COperator::Multi_MtV(const std::vector<CMatrix> &A,const std::vector<double> &b,std::vector<double> &vres)//9.SRes=(A^T)*b
{
	unsigned int i,j;

	unsigned int rowid,columnid;

	if(b.size()==vres.size())
	{
		for(i=0;i<vres.size();i++)
		{
			vres[i]=0.0;
		}

		for(i=0;i<A.size();i++)
		{
			columnid=A[i]._RowID;
			rowid=A[i]._ColID;

			vres[rowid-1]+=A[i]._Ele*b[columnid-1];
		}
	}
		
	else
	{
		cout<<"Error:结果向量的维数与相乘的矩阵和向量不同."<<endl;
	}
}

void COperator::Reve2ndMat(double a11,double a12,double a21,double a22,
								          double *ra11,double *ra12,double *ra21,double *ra22)//11.求二阶矩阵的逆矩阵
{
	double mod;

	mod=a11*a22-a12*a21;

	if(mod!=0)
	{
		*ra11= a22/mod; *ra12=-a12/mod;
		*ra21=-a21/mod; *ra22= a11/mod;
	}

	else
	{
		*ra11=1;
		*ra12=0;
		*ra21=0;
		*ra22=1;
	}
}

void COperator::Reve3rdMat(double a11,double a12,double a13,double a21,double a22,double a23,double a31,double a32,double a33,
													 double *ra11,double *ra12,double *ra13,double *ra21,double *ra22,double *ra23,double *ra31,double *ra32,double *ra33)//11.求三阶矩阵的逆矩阵
{
	double mod;

	mod=a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31);

	if(mod!=0)
	{
	 *ra11 = (a22*a33 - a23*a32)/mod;
	 *ra12 =-(a12*a33 - a13*a32)/mod;
	 *ra13 = (a12*a23 - a13*a22)/mod;
	 *ra21 =-(a21*a33 - a23*a31)/mod;  
	 *ra22 = (a11*a33 - a13*a31)/mod; 
	 *ra23 =-(a11*a23 - a13*a21)/mod;
	 *ra31 = (a21*a32 - a22*a31)/mod; 
	 *ra32 =-(a11*a32 - a12*a31)/mod;  
	 *ra33 = (a11*a22 - a12*a21)/mod;
	}
	else
	{
	 *ra11 = 1;
	 *ra12 = 0;
	 *ra13 = 0;
	 *ra21 = 0;  
	 *ra22 = 1; 
	 *ra23 = 0;
	 *ra31 = 0; 
	 *ra32 = 0;  
	 *ra33 = 1;
	}
}

double COperator::MinEig2ndMat(double a11,double a12,double a21,double a22)//求二阶矩阵的最小特征值
{
	double lambda1,lambda2;

	double temp1,temp2;

	temp1=a11*a11 - 2*a11*a22 + a22*a22 + 4*a12*a21;
	if(temp1<0)
		temp1=0.0;

	temp2=sqrt(temp1)/2;
 
	lambda1=a11/2 + a22/2 - temp2;
	 
	lambda2=a11/2 + a22/2 + temp2;

	if(lambda1>lambda2)
		return(lambda2);

	else
		return(lambda1);
}

double COperator::LMin( double a1,double a2 )//求两个数的最小值
{
	if (a1>a2)
	{
		return(a2);
	}

	else
		return(a1);
}

double COperator::LMax( double a1,double a2 )
{
	if (a1>a2)
	{
		return(a1);
	}

	else
		return(a2);
}




