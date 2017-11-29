#include"my_statistics.h"
#include<fstream>
#include<iostream>
#include"my_vrbls.h"
#include<numeric>
#include<math.h>
//Block->Data;
//Data will be analyzed in Analysis();
//------------Result-----------------
const int NObs = 5;
const int NQuan = 7;
std::vector<Data> Sample_Obs;
std::vector<Data> Block_Obs;
Data::Data(double ob0, double ob1, double ob2, double ob3, double ob4)
{
	Quan[0] = ob0;
	Quan[1] = ob1;
	Quan[2] = ob2;
	Quan[3] = ob3;
	Quan[4] = ob4;
}
/*
We calculate 0-2 in Nor_data();
0. The largest cluster size;
1. PercolatedQ;
2. Average Cluster Size;
3. S2
4. S4
5.Q1=<C1^2>/<C1>^2
6.Q2=<S2^2>/<3S2^2-2S4>
*/
char*  OutFile= "Percolation2D.dat";
char*  InFile = "Percolation2D_parameters.dat";
std::vector<double> Result_Ave;
std::vector<double> Result_Var;
std::vector<double> Result_Cor;
void Get_Parameters(std::istream &fin)
{
	std::cout << "Lx, Ly,Pc,NBlock,NSample,please:" << std::endl;
	fin >> Lx;
	fin >> Ly;
	fin >> Poc;
	fin >> NBlock;
	fin >> NSample;
	Block_Obs.resize(NBlock);
	Sample_Obs.resize(NSample);
}
void Get_Parameters(std::ifstream &fin)
{
	std::cout << "Lx, Ly,Pc,NBlock,NSample,please:" << std::endl;
	fin >> Lx;
	fin >> Ly;
	fin >> Poc;
	fin >> NBlock;
	fin >> NSample;
	Block_Obs.resize(NBlock);
	Sample_Obs.resize(NSample);
}
void write2file(std::ostream &fout)
{

	//fout << "Nblock	Nsample	Lx	Ly	p" << std::endl;
	fout << NBlock << "	" << NSample << "	" << Lx << "    " << Ly << "	 " << Poc<<" ";
	//fout << "		Result:		" << std::endl;
	//fout << "Lx	Ly	Poc	Quan	Ave	Var	Cor" << std::endl;
	for (int i = 0;i < NQuan;i++)
	{
		fout <<i << "	" << Result_Ave[i] << "	" << sqrt(Result_Var[i]) << "	" << Result_Cor[i] / Result_Var[i] <<" ";

	}
	fout<<std::endl;
}
void Collect_data(Block& p, int i)//collect date from block p into Sample_Obs(i)
{
	//i ranges from 0-NSample-1
	p.Calculate_Quan();
	for(int j=0;j<NObs;j++)
	{
	Sample_Obs[i].Quan[j] = p.Quan[j];
	}
}
void Normalize_data(int p) //Normalize date and write it into Block_Obs(p)
{
	double temp;
	double s2,c2;
	for (int j = 0;j < NObs;j++)
	{
		temp = 0;
		for (int i = 0;i < NSample;i++)
		{
			temp += Sample_Obs[i].Quan[j];
		}
		Block_Obs[p].Quan[j] = temp/NSample;
	}//calculate <C1>,<Q>,<S1>
	s2=0;c2=0;
	for(int i=0;i<NSample;i++)
	{
		c2+=Sample_Obs[i].Quan[0]*Sample_Obs[i].Quan[0];
		s2+=Sample_Obs[i].Quan[3]*Sample_Obs[i].Quan[3];
	}
	s2=s2/NSample;
	c2=c2/NSample;
	Block_Obs[p].Quan[5]=c2/(Block_Obs[p].Quan[0]*Block_Obs[p].Quan[0]);
	Block_Obs[p].Quan[6]=s2/(3*s2-2*Block_Obs[p].Quan[4]);

	/*
	We calculate 0-2 in Nor_data();
	0. The largest cluster size; (C1)
	1. PercolatedQ;
	2. Average Cluster Size;
	3. S2
	4. S4
	5.Q1=<C1^2>/<C1>^2
	6.Q2=<S2^2>/<3S2^2-2S4>
	*/
}
void Analyze_data()
{
	Result_Ave.resize(NQuan);
	Result_Var.resize(NQuan);
	Result_Cor.resize(NQuan);
	int i,j;
	double ave,var,cor;
	double devn, devp;  //deviation now and deviation prev
	for (i = 0;i < NQuan;i++)
	{
	//----------get ave-------------
		ave = 0;
		for (j = 0;j < NBlock;j++)
		{
			ave += Block_Obs[j].Quan[i];
		}
		ave = ave / NBlock;
		Result_Ave[i] = ave;
	//----------get ave finish-------------
	//----------get var and cor-------------
		devn = 0;devp = 0;cor = 0;var = 0;
		for (j = 0;j < NBlock;j++)
		{
			devn = Block_Obs[j].Quan[i] - ave;
			var += devn*devn;
			cor += devn*devp;
			devp = devn;
		}
		Result_Var[i] = var / NBlock;
		Result_Cor[i] = cor / NBlock;
	//----------get var and cor finish-------------
	}
}
