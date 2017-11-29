/*These are codes about random-number generator:
Most of it are copied from my_rng.f90.
If something makes wrong, you can check the original version
of this random-number generator.
Translated by Qiyang 2017.1.28
*/
#include <time.h>
#ifndef _MY_RNG
#define _MY_RNG
#ifndef _INC_MATH
#include <math.h>
#endif
#ifndef _STDIO_H_
#include <stdio.h>
#endif
/*
The initialization function of random-number generator
*/
int set_seed();
// for possible time-seed settings
extern int seed;
/*
rn() will return a double type random-number in U[0,1]
*/
double rn();
/*
rn_(i) will return an integer type random-number in [- 2147483647, 2147483647];
*/
int rn_i();
//-----------check elapsed time------------
extern time_t t_prev;
extern time_t t_now;
void set_elapse_time();
void elapse_time();

/**/
#endif
