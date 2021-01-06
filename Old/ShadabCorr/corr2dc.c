#include "corr2dc.h"

int cwrap(int x, int y) {

    // Cell/box wrap, required by periodic boundary conditions.
    //
    // Parameters:
    // x : the number to be wrapped.
    // y : the maximum/boundary of x.
    //
    // Returns:
    // modulo : the number insides the boundary.

    int modulo = x;
    if (x >= y) modulo = x - y;
    if (x < 0)  modulo = x + y;
    return modulo;
}

int relwrap(double *dxyz, double *blen) {

    // Relative wrap, required by periodic boundary conditions.
    //
    // Parameters:
    // dxyz : the number to be wrapped.
    // blen : the maximum/boundary of x.
    //
    // Returns:
    // modulo : the number insides the boundary.

    int ii;
    
    for(ii=0; ii<3; ii++) {
       if (dxyz[ii] >= blen[ii]/2) dxyz[ii] = dxyz[ii]-blen[ii];
       if (dxyz[ii] < -blen[ii]/2) dxyz[ii] = dxyz[ii]+blen[ii];
    }
    return 1;
}

double maxradius(double *rlim , int samp) {
   double maxrad;
   if(samp==0 || samp==2) { //rmu, rtheta
      maxrad=rlim[1];
   } else if (samp==1) { //rper_rpar
       double y;
       if(fabs(rlim[2])>fabs(rlim[3])) y=rlim[2];
       else y=rlim[3]; 
       maxrad=sqrt(rlim[1]*rlim[1]+y*y);
   } else if (samp==3) { //log(rper)-rpar
       maxrad=sqrt(pow(10,rlim[1]*2)+rlim[3]*rlim[3]);
   } else if (samp==4) { //log(r)-theta
       maxrad=pow(10,rlim[1]);
   } else {
      printf("\n****bad sampmode specified****\n");
      maxrad=-1;
      exit(0);
   }

   return maxrad;
}

void init_mesh(long *ll, long *hoc, double *p, long np, int nattr, int nho, double *blen,
      double *posmin){

    // Initialize the space. Mesh the particles in p into cells.
    //
    // Parameters:
    // ll : (np,) array, linked list of particles.
    // hoc : (nhocells, nhocells, nhocells) array, mesh of particles.
    // p : (np, nattr) array, particles.
    // np : number of particles.
    // nattr: number of features of particles.
    // nho : number of cells in one dimension.
    // blen : box length in 3d or the extent of survey.
    // posmin: lower extent of the survey

    long i;
    int ix, iy, iz, ii, jj, kk;

    // Initiate hoc to -1 for not missing the last point
    for (ii = 0; ii < nho; ii++)
        for (jj = 0; jj < nho; jj++)
            for (kk = 0; kk < nho; kk++)
                hoc[ii*nho*nho+jj*nho+kk] = -1;

    for (i = 0; i < np; i++) {
        ix = floor((p[i*nattr] - posmin[0]) / blen[0] * nho);
        iy = floor((p[i*nattr+1]-posmin[1]) / blen[1] * nho);
        iz = floor((p[i*nattr+2]-posmin[2]) / blen[2] * nho);

	if(ix==nho) ix=ix-1;
	if(iy==nho) iy=iy-1;
	if(iz==nho) iz=iz-1;

        ll[i] = hoc[ix*nho*nho+iy*nho+iz];
        hoc[ix*nho*nho+iy*nho+iz] = i;
    }
}

void corr2d(double *xc_2d, double *p1, long np1, double *p2, long np2, 
      double *rlim, int nbins0, int nbins1, int nhocells, double *blen, 
      double *posmin, int samp, int njn, int pbc, int los, int interactive) {

    // 2d correlation function.
    //
    // Parameters:
    // xc_2d : (nbins, nbins) array, 2d correlation to be returned.
    // p1 : (np1,6/7) first group of particles.
    // np1: number of particles in first group.
    // p2 : (np2,6/7) second group of particles.
    // np2 : number of particles of the second group.
    // rlim[4] : the minimum,maximum range you are looking at in both axis.
    // nbins0,1: number of bins in both axis.
    // nhocells: number of cells.
    // blen[3]: box size in 3d or the extent of data.
    // posmin[3]: minimum value of each co-ordinate (minimum corner of box)
    // samp : to decide the sampling scheme of the correlation function
    // njn : (0/njn) to decide if jacknife subsample needs to be evaluated or not
    // pc : (0/1) to decide periodic boundary condition or not

    long iq, i;
    int pp, qq, rr, p, q, r, ix, iy, iz; //To run loops on grids and particles
    int ppmin, ppmax, qqmin, qqmax, rrmin, rrmax; //To figure of the range of grid span
    int indx , indy; //To store the index of correlation function fbin
    double maxrad;
    int check, nattr;

    int nbins[2];
    nbins[0]=nbins0; nbins[1]=nbins1;
    //Allocate some memory
    //int *ncb = (int *)calloc(3, sizeof(int));
    //double *xyzp1 = (double *)calloc(3, sizeof(double));
    //double *xyzp2 = (double *)calloc(3, sizeof(double));
   
    int ncb[3], jn1,jn2,xcind;
    double xyzp1[3], xyzp2[3], weight12;
    long nvalid_pair=0, ninvalid_pair=0;

    maxrad=maxradius(rlim ,samp);

    for (pp=0; pp<3;pp++)
      ncb[pp] = ceil((maxrad / blen[pp]) * (double)(nhocells)) + 1;
      //ncb[pp] = floor((maxrad / blen[pp]) * (double)(nhocells)) + 1;

    long *ll = (long *)calloc(np2, sizeof(long));
    long *hoc = (long *)calloc(nhocells*nhocells*nhocells, sizeof(long));

    // Add weight or not
    if (njn > 0)  nattr =5 ;
    else  nattr = 4;

    // initialization
    init_mesh(ll, hoc, p2, np2, nattr, nhocells, blen, posmin);

    for (iq = 0; iq < np1; iq++) {
        if (iq % 100000 == 0 && interactive>=0) printf(".");

        //read in first particle co-ordinate
        xyzp1[0] = p1[iq*nattr]; xyzp1[1] = p1[iq*nattr+1]; xyzp1[2] = p1[iq*nattr+2];

        ix = floor((xyzp1[0]-posmin[0]) / blen[0] * nhocells);
        iy = floor((xyzp1[1]-posmin[1]) / blen[1] * nhocells);
        iz = floor((xyzp1[2]-posmin[2]) / blen[2] * nhocells);

        
	ppmin=ix-ncb[0]; ppmax=ix+ncb[0];
	qqmin=iy-ncb[1]; qqmax=iy+ncb[1];
	rrmin=iz-ncb[2]; rrmax=iz+ncb[2];
	
	if(pbc!=1) { //pbc can be applied only with los=1 along z axis
           if(ppmin<0) ppmin=0;
           if(qqmin<0) qqmin=0;
           if(rrmin<0) rrmin=0;

           if(ppmax>nhocells) ppmax=nhocells;
           if(qqmax>nhocells) qqmax=nhocells;
           if(rrmax>nhocells) rrmax=nhocells;
	}

        for (pp = ppmin; pp < ppmax; pp++) {
            p = cwrap(pp, nhocells);
            for (qq = qqmin; qq < qqmax; qq++) {
                q = cwrap(qq, nhocells);
                for (rr = rrmin; rr < rrmax; rr++) {
                    r = cwrap(rr, nhocells);
                    if (hoc[p*nhocells*nhocells+q*nhocells+r] != -1) {
                        i = hoc[p*nhocells*nhocells+q*nhocells+r];
                        while (1) { //while follow the chains
                            //read in second particle co-ordinate
                            xyzp2[0]=p2[i*nattr]; 
               			    xyzp2[1]=p2[i*nattr+1]; 
			                   xyzp2[2]=p2[i*nattr+2];
                            weight12=p1[iq*nattr+3]*p2[i*nattr+3];
			     
                            check=find_bin(xyzp1, xyzp2,samp,nbins,rlim, 
				                los,pbc,blen,&indx,&indy);

                            if(check!=-1 && 
				                 indx< nbins[0] && indy < nbins[1] &&
				                 indx>=0 && indy>=0) {
			                      if(njn==0) { //Taking care of jacknife regions
				                       xcind=indx*nbins[1]+indy;
			                          xc_2d[xcind] = xc_2d[xcind]+weight12;
			                      }
			                      else { //else-jacknife
				                       xcind=indx*nbins[1]*(njn+1)+indy*(njn+1)+njn;
			                          xc_2d[xcind] = xc_2d[xcind]+weight12;
				                       //Jacknife region of particle 1
				                       jn1=p1[iq*nattr+4]; 
				                       xcind=indx*nbins[1]*(njn+1)+indy*(njn+1)+jn1;
			                          xc_2d[xcind] = xc_2d[xcind]+weight12;
				                       //Jacknife region of particle 2
				                       jn2=p2[i*nattr+4]; 
                                   xcind=indx*nbins[1]*(njn+1)+indy*(njn+1)+jn2;
				                       if(jn1!=jn2) xc_2d[xcind] = xc_2d[xcind]+weight12;
			                      } //end else-jacknife
			                      nvalid_pair++;
			                   } //end if-check
			                   else
			                      ninvalid_pair++;

                            if (ll[i] != -1) {
                               i = ll[i];
                            }
                            else break; //This breaks the while
                        } //This is end of while(1)
                    } //This is the end of hoc conditions
                } //This is the end of rr for-loop
            } //This is the end of qq for-loop
        }//This is the end of pp for-loop
    } //This is the end of iq for-loop
    free(ll);
    free(hoc);
    if(interactive>1)
       printf("\nnvalid_pair= %ld , invalid_pair= %ld \n",nvalid_pair, ninvalid_pair);
} //This is the end of function corr2d

//function to define the sampling and the bin
int find_bin(double *axyz, double *bxyz, int samp, int *nsamp, double *sampbound, 
      int los, int pbc, double *blen,int *rx, int *ry) {
   double difxyz[3], difmag=0, difDOTlos=0;
   double losxyz[3], losmag=0;
   double xsep=0,ysep=0, xcoord=0, ycoord=0;
   int ii;

   //Defining the line of sight
   if (los==0 || los==2 || los==3) {// the los is along the mid/first/second galaxy co-ordinate
      for (ii=0;ii<3;ii++){
         difxyz[ii]=axyz[ii]-bxyz[ii];
	 if(los==0) losxyz[ii]=0.5*(axyz[ii]+bxyz[ii]); //los is along the mid point
	 else if(los==2) losxyz[ii]=axyz[ii];//los is along the first galaxy
	 else losxyz[ii]=bxyz[ii]; // los is along the second galaxy
         losmag=losmag+losxyz[ii]*losxyz[ii];
         difmag=difmag+difxyz[ii]*difxyz[ii];
         difDOTlos=difDOTlos+difxyz[ii]*losxyz[ii];
      }
      difmag=sqrt(difmag);
      losmag=sqrt(losmag);
      ycoord=difDOTlos/losmag;  
      xcoord=sqrt(difmag*difmag-ycoord*ycoord);
   }
   else if(los==1) {//the los is along the Z axis
      for (ii=0;ii<3;ii++)
         difxyz[ii]=axyz[ii]-bxyz[ii];
      if(pbc==1) { //pbc can be applied only with los=1 along z axis
	 relwrap(difxyz, blen); 
      }
      ycoord=difxyz[2];
      xcoord=sqrt(difxyz[0]*difxyz[0]+difxyz[1]*difxyz[1]);
   }

   //Using a given sampling definition for 2d correlation function
   if(samp==0) { //sampling in rmu space
      xsep=sqrt(xcoord*xcoord+ycoord*ycoord);
      ysep=fabs(ycoord/xsep); //consider positive and negative mu together
   }
   else if (samp==1) { //rpar-rper
      ysep=ycoord;
      xsep=xcoord;
   }
   else if (samp==2) { //r theta
      xsep=sqrt(xcoord*xcoord+ycoord*ycoord);
      ysep=acos(ycoord/xsep);
   }
   else if (samp==3) { //log rperp- rpar
      ysep=ycoord;
      xsep=log10(xcoord);
   }
   else if (samp==4) { //log(r)-theta
      xsep=sqrt(xcoord*xcoord+ycoord*ycoord);
      ysep=acos(ycoord/xsep);
      xsep=log10(xsep);
   }
   else {
      printf("\ninvalid sampling mode %d\n",samp);
      exit(0);
   }

   //Determine the bin in which the pair lies
   if(xsep > sampbound[0] && xsep < sampbound[1] && ysep >sampbound[2] && ysep <sampbound[3])
   {
      int xbin_n = floor((xsep-sampbound[0])*nsamp[0]/(sampbound[1]-sampbound[0]));
      int ybin_n = floor((ysep-sampbound[2])*nsamp[1]/(sampbound[3]-sampbound[2]));
      if(xbin_n == nsamp[0]) xbin_n -= 1;
      if(ybin_n == nsamp[1]) ybin_n -= 1;
      *rx=xbin_n;
      *ry=ybin_n;
      return 0;
   }
   else {
      *rx=-1;
      *ry=-1;
      return -1;
   }

   return -1;
}

