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

int main (int argc, char **argv) {

  FILE * DD;
  FILE * PP;
  FILE * RR;

  const int Nr = 150;
  const int Nm = 100;
  const double RMIN = 0.0;
  const double RMAX = 150.0;
  const double RSTEP = 1.0;
  const double MSTEP = 0.01;

  const int maxNgal = 40000000;

  int dummy;

  char DDname[400];
  char PPname[400];
  char RRname[400];
  
  int i;
  int bit;
  int bittot;

  int *n;

  sprintf(DDname,argv[1]);
  sprintf(RRname,argv[2]);
  sprintf(PPname,argv[3]);

  struct galaxy * gald;
  struct galaxy * galr;

  double * DRcount;
  double * DRcount_in_place;
  gald = (struct galaxy*)malloc(maxNgal*sizeof(struct galaxy));
  galr = (struct galaxy*)malloc(maxNgal*sizeof(struct galaxy));
  DRcount = (double*)calloc(Nr*Nm,sizeof(double));
  DRcount_in_place = (double*)calloc(Nr*Nm,sizeof(double));
  struct galaxy * galpd, * galpr;
  galpd = gald;
  galpr = galr;

  DD = fopen(DDname,"r");
  RR = fopen(RRname,"r");

  double ra, dec, red, dist;
  double dd;
  int Nd = 0;
  double weight;
  double ddummy;
  int indexFKP;
  double weightFKP;
  double weight1, weight2, weight3, weight4;
  char line [2000];

  int N = 0;
  double sumW = 0;
  double sumW_in_place = 0;

  for (int i = 0; i < 0; i ++) {
  fgets (line, 2000, DD);
  printf("%s\n", line);
  }
  while (fscanf(DD,"%lf %lf %lf\n",&ra,&dec,&red) != EOF) {
    ra *= PI/180.0;
    dec *= PI/180.0;
    dist = rz(red);
    galpd->x = dist*cos(dec)*cos(ra);
    galpd->y = dist*cos(dec)*sin(ra);
    galpd->z = dist*sin(dec);
    galpd->d = dist;
//    galpd->w = weight1 * weight2;
    galpd++;
    Nd++;
  }
  fclose(DD);

  N = 0;
  for (int i = 0; i < 0; i ++) {
  fgets (line, 2000, RR);
  printf("%s\n", line);
  }
  while (fscanf(RR,"%le %le %le\n",&ra,&dec,&red) != EOF) {
    ra *= PI/180.0;
    dec *= PI/180.0;
    dist = rz(red);
    galpr->x = dist*cos(dec)*cos(ra);
    galpr->y = dist*cos(dec)*sin(ra);
    galpr->z = dist*sin(dec);
    galpr->d = dist;
//    galpr->w = weight1 * weight2;
    galpr++;
    N++;
  }
  fclose(RR);  

  MPI_Init (&argc, &argv);
  int numproc;
  int totproc;
  MPI_Comm_rank (MPI_COMM_WORLD, &numproc);
  MPI_Comm_size (MPI_COMM_WORLD, &totproc);
  
  n = (int *) calloc(totproc + 1, sizeof(int));
  n[0] = 1;
  for (i = 1; i < totproc; i++) {
    n[i] = (int)floor(double(N)/totproc*i);
  }
  n[totproc] = N;
  bit = n[numproc];
  bittot = n[numproc + 1];
  

  sprintf(PPname, "%s.DR",PPname);
  double mu;
  double d1, d2, d3;
  double x1, y1, z1, x2, y2, z2, gd1, gd2, w1, w2;
  double r2, rat;
  double xb, yb, zb, db2;
  double rr;
  galpd = gald;
  int j;

  for (i = 0; i < Nd; i++) {
    x1 = galpd->x; y1 = galpd->y; z1 = galpd->z; gd1 = galpd->d; //w1 = galpd->w;
    galpr = galr + bit;
    for (j = bit; j < bittot; j++) {
       x2 = galpr->x; y2 = galpr->y; z2 = galpr->z; gd2 = galpr->d; //w2 = galpr->w;
       //sumW += w1*w2;
       sumW += 1.0;
       d1 = x1 - x2;
       d2 = y1 - y2;
       d3 = z1 - z2;
       r2 = d1*d1 + d2*d2 + d3*d3;
       if (r2 > RMAX*RMAX || r2 < RMIN*RMIN) {
	 galpr++;
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
//         DRcount[ind] += w1*w2;
         DRcount[ind] += 1.0;
       }
       galpr++;
    }
    galpd++;
  }
  free(gald);
  free(galr);
  free(n);

  MPI_Reduce(DRcount, DRcount_in_place, Nr*Nm, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumW, &sumW_in_place, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  printf("weighted number = %le %d\n", sumW, numproc);
  if (numproc == 0) {
    PP = fopen(PPname,"w");
    printf("%s\n", PPname);
    fprintf(PP,"# weighted number: %lf\n",sumW_in_place);
    fprintf(PP,"# RBINS: %d\n", Nr);
    fprintf(PP,"# MBINS: %d\n", Nm);
    int k, l;
    for (k = 0; k < Nm; k++) {
      for (l = 0; l < Nr; l++) {
      fprintf(PP,"%lf ",DRcount_in_place[k*Nr + l]);
      }
      fprintf(PP,"\n");
    }
    fclose(PP);
  }

  free(DRcount);
  free(DRcount_in_place);

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
