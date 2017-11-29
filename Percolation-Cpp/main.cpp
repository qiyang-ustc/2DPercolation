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
	Vol = Lx*Ly;
	Block_Obs.resize(NBlock);
	Sample_Obs.resize(NSample);
	TotalSample = NBlock*NSample;
	label.Resize(Vol);
}
void GenerateBlock(Block b)
{
	int i;
	int j;
	int left;
	int up;
	ident = 1;
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

				if (up == 0 && left == 0)
				{
					b.ele(i, j) = ident;
					ident++;
				}
				else if (up == 0 && left!=0)
				{
					b.ele(i, j) =label.getfather(up+left);
				}
				else if (left == 0 && up != 0)
				{
					b.ele(i, j) = label.getfather(up + left);
				}
				else
				{
					b.ele(i, j) = label.getfather(left);
					label.Union(up, left);
				}

			}
			else
			{
				b.ele(i, j) = 0;
			}
		}

	}
	for(i=0;i<Lx;i++)
		for (j = 0;j < Ly;j++)
		{
			if (b.ele(i, j) != 0)
			{
				b.ele(i, j) = label.getfather(b.ele(i, j));
			}
		}
	label.Refresh();
}
time_t tm;

int main(int argc,char *argv[])
{
	Lx = atoi(argv[1]);
	Ly = atoi(argv[2]);
	Poc=atof(argv[3]);
	NBlock=atoi(argv[4]);
	NSample=atoi(argv[5]);
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
	set_elapse_time();
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
