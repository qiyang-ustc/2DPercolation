
//*******************************************************************
// Ising model on the square Lattice

// Error bars are calculated using the blocking technique.
// Let ' \chi^2 = \sum_{t=1}^{T} <(O_t -<O>)^2> ' be the squared chi for
// 'T' blocks of observable 'O'. Assuming each block of data is independent
// of each other, the error bar of 'O' is given by 'Dev = \sqrt{\chi^2/T(T-1)}.

// Reliabity of the obtained errors is monitored by t=1 correlation,
// for which tolerance is set by variable 'tol' (default: tol=0.20d0).

// Composite quantities like Binder ratios are calculated in each block, and
// the associated error bars are obtained from their fluctuations.

// Results are written into a special file 'dat.***' if the number of
// blocks is less than 125 or correlation is too big. Data in each
// block will be also printed out in this case.

// Default number of extensive simulation is 'NBlck=1024'.

// For test purpose, for which huge amount of information will be
// printed out, 'NBlck' should be set smaller but >2.

// Dynamical behavior is not studied.

// 'my_vrbls.f90', 'carlo.f90', 'monte.f90', 'measure.f90',
// 'write2file.f90' and etc need to be modified for new projects.

//  Look for 'PROJECT-DEPENDENT'.

//  Author: Yuan Huang
//  Date  : April 19th, 2012.
//*******************************************************************

// Look for 'PROJECT-DEPENDENT' for different projects
#include <vector>
#include<iostream>
#include"my_rng.h"
#ifndef _MY_VRBLS_
#define _MY_VRBLS_
//------------------Project-Independed---------
extern const double tm32;	 //1/2^(32)
extern const double eps;     //very small number
extern const double tol;	 //tolerance for correlation 0.2
extern const int MaxInt;	 //Maximun integer
extern const int MinInt;	 //Minimun integer
extern const int MaxBlock;	     //Maximun number of block
extern const int MinBlock;		 //Minimun number of block
extern int NBlock;               //N=number of
extern int NSample;				// samples in a block
extern int TotalSample;
extern int NToss;				//samples to be trown away
extern int Collect_Interval;
//----------------Project-depended-------------
extern const int Dimension;
extern int Lx, Ly;
extern double Poc;
extern int Vol;
extern int ident;//used to identify cluster
extern const int NObs;
extern const int NQuan;
template<class T>
/* Some algorithm which may be useful*/
void Delete_Num(typename std::vector<T>& v, T del)
//Can not Delete Points!!!
{
	for (typename std::vector<T>::iterator iter = v.begin();iter != v.end();)
	{
		if (del == *iter)
		{
			iter = v.erase(iter);
		}
		else
		{
			++iter;
		}
	}
}
//----------ALGORITHM------END!!!

//---------------Block-----------------------
class Block {
public:
	std::vector<double> Quan;
	Block()
	{
		this->p_block = new int*[Lx];
		for (int i = 0;i < Lx;i++)
		{
			(this->p_block)[i] = new int[Ly];
		}
		temp.resize(Vol);
		for (int i = 0;i < Vol;i++)
		{
			temp[i] = 0;
		}
		percolated_num.resize(Vol);
		Quan.resize(NObs);
		//initialize temp  which is the hash function.
	}
	~Block() {}
	void Print_ele(std::ostream &fout)
	{
		for (int i = 0;i < Lx;i++)
		{
			for (int j = 0;j < Ly;j++)
			{
				fout << this->ele(i, j) << "	";
			}
			fout << std::endl;
		}
	}//out put block configuration.
	int& ele(int i, int j)
	{
		return (this->p_block)[i][j];
	}//return p_block(i,j)
	void replace_a2b_UNTILcd(int a, int b,int c,int d)
	{//substitude a to b;
		for(int i=0;i<=c;i++)
			for (int j = 0;j < Ly;j++)
			{
				if ((this->p_block)[i][j] == a)
					this->p_block[i][j] = b;
			}
	}
	void fresh()
	{
		for (int i = 0;i < Lx;i++)
			for (int j = 0;j < Ly;j++)
				p_block[i][j] = 0;
	}
	int Occupied_Number()
	{
		int num=0;
		for (int i = 0;i < Lx;i++)
		{
			for (int j = 0;j < Ly;j++)
			{
				if (this->ele(i, j) != 0)
				{
					num++;
				}
			}
		}
		return num;
	}
	int PercolatedQ()
	{
		std::vector<int> up,down;
		int Q=0;
		bool temp;
		for (int i = 0;i < Vol;i++)
		{
			percolated_num[i] = 0;
		}
		for (int i = 0;i < Ly;i++)
		{
			if (this->ele(0, i) != 0)
				up.push_back(this->ele(0, i));
			if (this->ele(Ly-1, i) != 0)
				down.push_back(this->ele(Ly-1, i));
		}
		while(up.size()!=0)
		{
			temp=0;
			for(int i=0;i<down.size();i++)
			{
				if(up.front()==down[i])
				{temp=1;
				break;}
			}
			if (temp)
			{
				Q = 1;
				percolated_num[*(up.begin())]++;
				Delete_Num(up, *(up.begin()));
			}
			else
			{
				Delete_Num(up, *(up.begin()));
			}
		}
		//this->Print_ele(std::cout);
		//std::cout << Q;
		if (Q == 1)
			return 1;
		else
			return 0;
	}
	double Average_Cluster_Size()
	{
		long c1=0;long s2=0;long s4=0;
		int cluster_num = 0;
		int site_num = 0;
		for (int i = 0;i < Lx;i++)
		{
			for (int j = 0;j < Ly;j++)
			{
				if (this->ele(i, j) != 0)
				{
					temp[ele(i,j)]++;
				}
			}
		}

		//-----refresh block and we need to calculate cluster number.
		for (int i = 0;i < Vol;i++)
		{
			if (temp[i] != 0)
			{
				if(percolated_num[i]==0)
				{
					site_num += temp[i];
					cluster_num++;
				}
			if(temp[i]>c1)
			{
				c1=temp[i];
			}
			s2+=temp[i]*temp[i];
			s4+=s2*s2;
				temp[i] = 0;
			}
		}
		this->Quan[0]=c1;
		this->Quan[3]=s2;
		this->Quan[4]=s4;
		this->Occupied_Sites = site_num;
		return (double)site_num / cluster_num;
	}
	void Calculate_Quan()
	{
		/*
		We calculate 0-2 in Nor_data();
		0. The largest cluster size; (C1)
		1. PercolatedQ;
		2. Average Cluster Size;
		3. S2
		4. S4
		*/
		this->Quan[1]=this->PercolatedQ();
		this->Quan[2]=this->Average_Cluster_Size();
		//calculation of 0,3,4  are finished in Average_Cluster_Size();

	}
private:
	int **p_block;
	std::vector<int> temp;
	int Occupied_Sites;
	std::vector<int> percolated_num;

};




#endif
