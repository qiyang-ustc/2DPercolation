#include <vector>
bool intersected(std::vector<int> up,std::vector<int> down)// return 1 if a ann b are intersected.
{
	int temp=0;
	for(int i=0;i<up.size();i++)
		for(int j=0;j<down.size();j++)
		{
				if(up[i]==down[j])
					return 1;
		}
   return 0;
}
