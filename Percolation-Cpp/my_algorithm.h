#include <vector>
#include<iostream>
#include"my_rng.h"
#ifndef _MY_ALGO_
#define _MY_ALGO_
class Union_Find_int
{
public:
	Union_Find_int(int size)
	{
		labels.resize(size);
		for (int i = 0;i < size;i++)
		{
			labels[i] = i;
		}
	}
	Union_Find_int() {}
	void Refresh() 
	{
		int size = labels.size();
		for (int i = 0;i < size;i++)
		{
			labels[i] = i;
		}
	}
	int getfather(int x)
	{
		int y=x;
		int z;
		while (y != labels[y])
		{
			y = labels[y];
		}
		while (x != labels[x])
		{
			z = labels[x];
			labels[x] = y;
			x = z;
		}//compress path;
		return y;
	}
	void Union(int x, int y)
	{
		if(x!=y)
		labels[getfather(x)] = getfather(y);
	}
	void Resize(int t)
	{
		labels.resize(t);
		this->Refresh();
	}
private:
	std::vector<int> labels;
};

#endif
