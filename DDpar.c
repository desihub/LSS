#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_histogram.h>
#include"mpi.h"

const double PI = 3.14159;

double rz(double red);

struct galaxy {
  double x, y, z, d, w;
};

void N2divide( int * n, int N, int Ntot);

int main (int argc, char **argv) {


  printf("starting program\n");
  FILE * DD;
  FILE * PP;

  const int Nr = 150;
  const int Nm = 100;
  const double RSTEP = 1.0;
  const double MSTEP = 0.01;

  const int maxNgal = 3000000;
  const double RMIN = 0.0;
  const double RMAX = 150.0;

  int *n;

  char DDname[400];
  char PPname[400];
  
  int bit;
  int bittot;

  sprintf(DDname,argv[1]);

  MPI_Init (&argc, &argv);
  int numproc;
  int totproc;
  MPI_Comm_rank (MPI_COMM_WORLD, &numproc);
  MPI_Comm_size (MPI_COMM_WORLD, &totproc);

  n = (int *) calloc(totproc + 1, sizeof(int));

  sprintf(PPname,"%s.DD",argv[2]);

  DD = fopen(DDname,"r");
  printf("%s %d\n",DDname,numproc);
  struct galaxy * gal;
  double * DDcount;
  double * DDcount_in_place;
  gal = (struct galaxy*)malloc(maxNgal*sizeof(struct galaxy));
  DDcount = (double*)calloc(Nr*Nm,sizeof(double));
  DDcount_in_place = (double*)calloc(Nr*Nm,sizeof(double));
  struct galaxy * galp;
  galp = gal;

  int N = 0;
  int Ntest = 0;
  double sumW = 0;
  double sumW_in_place = 0;
  double ra, dec, red, dist;
  double weight = 1.0;
  double weight1;
  double weight2;
  double weight3;
  double weight4;
  double weightFKP;
  int indexFKP;
  double ddummy;
  char line[2000];
  for (int i = 0; i < 0; i ++) {
  fgets (line, 2000, DD);
  printf("%s\n", line);
  }
  printf("start reading catalogue\n");
  while (fscanf(DD,"%lf %lf %lf\n",&ra,&dec,&red) != EOF) {
    ra *= PI/180.0;
    dec *= PI/180.0;
    dist = rz(red);
    galp->x = dist*cos(dec)*cos(ra);
    galp->y = dist*cos(dec)*sin(ra);
    galp->z = dist*sin(dec);
    galp->d = dist;
   // galp->w = weight1 * weight2;
    galp++;
    N++;
  }
  fclose(DD);
  n[0] = 1;
  N2divide(n, totproc, N);
  n[totproc] = N;
  n[0] = 0;
  bit = n[numproc];
  bittot = n[numproc + 1];
  for (int i = 0; i < totproc + 1; i++) {
    printf("n[%d] = %d\n", i, n[i]);
  } 
  double mu;
  double d1, d2, d3;
  double x1, y1, z1, x2, y2, z2, gd1, gd2;// gw1, gw2;
  double r2, rat;
  double rr;
  double xb, yb, zb, db2;
  struct galaxy * galp1, * galp2;  
  galp1 = gal;
  galp1 += bit;
  for (int i = bit; i < bittot; i++) {
    x1 = galp1->x; y1 = galp1->y; z1 = galp1->z; gd1 = galp1->d; //gw1 = galp1->w;
    galp2 = gal + i + 1;
    for (int j = i + 1; j < N; j++) {
       x2 = galp2->x; y2 = galp2->y; z2 = galp2->z; gd2 = galp2->d; //gw2 = galp2->w;
    //   sumW += gw1*gw2;
       sumW += 1.0;
       d1 = x1 - x2;
       d2 = y1 - y2;
       d3 = z1 - z2;
       r2 = d1*d1 + d2*d2 + d3*d3;
       if (r2 > RMAX*RMAX || r2 < RMIN*RMIN) {
	   galp2++;
         continue;
       }
       rat = gd1/gd2;
       xb = x1 + x2*rat;
       yb = y1 + y2*rat;
       zb = z1 + z2*rat;
       db2 = xb*xb + yb*yb + zb*zb;
       mu = fabs((xb*d1 + yb*d2 + zb*d3)/sqrt(r2)/sqrt(db2));
       rr = sqrt(r2);
	int binr = (int)((rr - RMIN)/RSTEP);
	int binm = (int)(mu/MSTEP);
	if (binr >= 0 && binm >= 0 && binr < Nr && binm < Nm){
	int ind = binr + Nr*binm;
//	DDcount[ind] += gw1*gw2;
        DDcount[ind] += 1.0;
	}
       galp2++;
    }
    galp1++;
  }
  free(gal);

  MPI_Reduce(DDcount, DDcount_in_place, Nr*Nm, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumW, &sumW_in_place, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (numproc == 0) {
    PP = fopen(PPname,"w");
    printf("%s\n", PPname);
    fprintf(PP,"# weighted number: %lf\n",sumW_in_place);
    fprintf(PP,"# RBINS: %d\n",Nr);
    fprintf(PP,"# MBINS: %d\n",Nm);
    int k, l;
    for (k = 0; k < Nm; k++) {
      for (l = 0; l < Nr; l++) {
      fprintf(PP,"%lf ",DDcount_in_place[k*Nr + l]);
      }
      fprintf(PP,"\n");
    }
    fclose(PP);
  }

  free(DDcount);
  free(DDcount_in_place);

  free(n);
  MPI_Finalize ();

  return 0;
}

double f (double x, void * p) {
  double Om = 0.310;
  double ff = 2997.92458/sqrt(Om*(1.0 + x)*(1.0 + x)*(1.0 + x) + 1.0 - Om);
  return ff;
}

double rz (double red) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  double result, error;
  gsl_function F;
  F.function = &f;
  gsl_integration_qags(&F, 0, red, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result;
}


void N2divide( int * n, int N, int Ntot) {
   double dummy;
   for (int i = 1; i < N; i++) {
     dummy = double(Ntot) - 0.5;
     dummy -=
sqrt((2*double(Ntot)-1)*(2*double(Ntot)-1)-4*((2*double(Ntot)-n[i-1])*(n[i-1]-1)+double(Ntot)*(double(Ntot)-1)/N))/2.0;
     n[i] = int(dummy);
   }

   return;

}
