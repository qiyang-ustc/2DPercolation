#include <iostream>
#include <fstream>
#include"my_vrbls.h"
#include<string.h>
#ifndef  __MY_STA__
#define __MY_STA__

//-----------Class Data------------
class Data
{
public:
	Data() {}
	Data(double, double, double,double,double);
	~Data() {};
	double Quan[7];
};
//----------------------------------
extern char* OutFile;
extern char* InFile;
extern void Get_Parameters(std::ifstream&);
extern void Get_Parameters(std::istream&);
extern void write2file(std::ostream&);
extern std::vector<Data> Sample_Obs;
extern std::vector<Data> Block_Obs;
extern void Collect_data(Block&, int);
extern void Normalize_data(int);
extern void Analyze_data();
#endif // ! __MY_IO__
