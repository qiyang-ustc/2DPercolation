# 2DPercolation

Only for site percolation.

There are 7 quantities which are measured in this program.
	0. The largest cluster size C1
	1. PercolatedQ;(Probablity)
	2. Average Cluster Size;
	3. S2:=\Sum size^2
	4. S4:=\Sum size^2
	5.Q1=<C1^2>/<C1>^2
	6.Q2=<S2^2>/<3S2^2-2S4>


Input:
Lx Ly Poc NBlock NSample
for example in Linux System:
./a.out 128 128 0.55 1024 1024

Output:
This program will put out the result into ./Percolation.dat
in order:
NBlock NSample Lx Ly poc 
and quantities:
Quantity_Number  Average  Standrad_Deviation	Correlation(useless in this model)
for example
100	100	256    256	 0.592746 0	15909.7	549.245	-0.0400003 1	0.4913	0.0476792	-0.0677677 2	15.3108	0.507715	-0.0145927 3	3.26786e+08	1.63047e+07	-0.0301495 4	1.90752e+18	3.83134e+17	0.0799993 5	1.13033	0.0154924	-0.204998 6	-0.0420523	0.0122523	0.029181 
