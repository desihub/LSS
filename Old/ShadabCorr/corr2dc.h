#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int cwrap(int x, int y);
int relwrap(double *dxyz, double *blen);
double maxradius(double *rlim , int samp);
void init_mesh(long *ll, long *hoc, double *p, long np, int nattr, int nho, double *blen, double *posmin);
void corr2d(double *xc_2d, double *p1, long np1, double *p2, long np2, double *rlim,
        int nbins0, int nbins1, int nhocells, double *blen, double *posmin, int samp, 
	int njn, int pbc, int los, int interactive);
int find_bin(double *axyz, double *bxyz, int samp, int *nsamp, 
      double *sampbound, int los, int pbc, double *blen, int *rx, int *ry);
