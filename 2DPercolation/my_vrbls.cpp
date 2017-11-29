#include <math.h>
#include <vector>
#include "my_rng.h"
#include <fstream>
const int      Mxint = 2147483647;			// maximum int
const int      Mnint = -2147483647;	    // minimum int
const double   tm32 = 1.0 / pow(2, 32);
const double   eps = 1.0E-14;				//very small number
const double   tol = 0.20;					//tolerance of correlation
const int      MaxInt = 2147483647;
const int      MinInt = -2147483647;
int NBlock;
int NSample;
int TotalSample;
int NToss;
int Collect_Interval;
int TossThrown;



int ident = 1;
const int Dimension = 2;
int Lx, Ly, Lz;
int Vol;
double Poc;
//------------Class spin and block--------------


