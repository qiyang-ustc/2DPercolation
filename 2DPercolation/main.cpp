#include<iostream>
#include<fstream>
#include"my_rng.h"
#include"my_vrbls.h"
#include"my_statistics.h"
#include"my_algorithm.h"
#include<vector>
#include<stdlib.h>
Union_Find_int label;
void Initialize()
{
	Vol = Lx*Ly;Vol2 = Vol*Vol;Vol4 = Vol2*Vol2;
	Block_Obs.resize(NBlock);
	Sample_Obs.resize(NSample);
	TotalSample = NBlock*NSample;
	label.Resize(Vol);
}
void GenerateBlock(Block b)/*for site and bond mixed percolation*/
{
	int i;
	int j;
	int left;
	int up;
	ident = 1;
	bool bondup;
	bool bondleft;
	for (i = 0;i < Lx;i++)
	{
		for (j = 0;j < Ly;j++)
		{
			if (rn() < Poc)
			{
				if (i == 0) { up = 0; }
				else { up = b.ele(i - 1, j); }
				if (j == 0) { left = 0; }
				else { left = b.ele(i, j - 1); }
				//generate bond
				bondup = (rn() < Pbc&&up != 0);
				bondleft = (rn() < Pbc&&left != 0);

				if (bondup == 1 && bondleft == 0)
				{
					b.ele(i, j) = label.getfather(up);
				}
				else if (bondup == 0 && bondleft == 1)
				{
					b.ele(i, j) = label.getfather(left);
				}
				else if (bondup == 1 && bondleft == 1)
				{
					b.ele(i, j) = label.getfather(left);
					label.Union(up, left);
				}
				else
				{
					b.ele(i, j) = ident;
					ident++;
				}
			}
			else
			{
				b.ele(i, j) = 0;
			}
		}

	}
//Now apply period boundary condition on it;
	/*
for(int i=0;i<Lx;i++)
{
	if(rn()<Pbc && b.ele(0,i)!=b.ele(Ly-1,i) && b.ele(Ly-1,i)!=0 && b.ele(0,i)!=0)
	{
		label.Union(b.ele(0,i),b.ele(Ly-1,i));
	}
}
for(int i=0;i<Ly;i++)
{
	if(rn()<Pbc && b.ele(i,0)!=b.ele(i,Lx-1)  && b.ele(i,0)!=0 && b.ele(i,Lx-1)!=0)
	{
		label.Union(b.ele(i,0),b.ele(i,Lx-1));
	}
}
	*/

//applying end
	for(i=0;i<Lx;i++)
	{for (j = 0;j < Ly;j++)
		{
			if (b.ele(i, j) != 0)
			{
				b.ele(i, j) = label.getfather(b.ele(i, j));
			}
		}
  }
	label.Refresh();
}
time_t tm;

int main(int argc,char *argv[])
{
	set_elapse_time();
	Lx = atoi(argv[1]);
	Ly = atoi(argv[2]);
	Poc=atof(argv[3]);
	Pbc=atof(argv[4]);
	NBlock=atoi(argv[5]);
	NSample=atoi(argv[6]);
	std::ofstream fout(OutFile, std::ios::app);
	//We use System parameters to get parameters;
	/*
	std::ifstream fin(InFile, std::ios::in);
	if (fin.is_open())
	{
		Get_Parameters(fin);
	}
	else
	{
		Get_Parameters(std::cin);
	}
	*/
	Initialize();
	Block b;
	for (int i = 0;i < NBlock;i++)
	{
		for (int j = 0;j < NSample;j++)
		{
			b.fresh();      //As we haven't fresh it in Average_cluster_size;
			GenerateBlock(b);
			Collect_data(b,j);
		}
		Normalize_data(i);
	}
	//-------Analysis----------
	Analyze_data();
	//-------write to file-----
	write2file(fout);
	elapse_time();
	//system("pause");
}
