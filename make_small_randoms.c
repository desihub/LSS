/** 
 * This programm will take a mask file - argv[1], an input catalogue - argv[2],
 * and will generate N = argv[4] number of random catalogues, with names -
 * argv[3] + xx where xx is a number from 00 to N, for a target of type -
 * argv[5].
 * 
 * The targets for which the random is made will have a random angular position
 * drawn from the mask file and a random redshift drawn from the target
 * catalogue (for the same type of objects). Targets of other type will be in
 * exactly the same position as in the original target catalogue.
 *
 * The random catalogue will have the same size as the input catalogue. The
 * random files will be run through a fiber-assignment code and concatenated.
 * The random file will have an object ID.
 *
 * The random catalogue will be binary (this is what Martin/Bob's current fiber
 * assignment code takes in. With the format: (int)Nobj, (float)ra[Nobj),
 * (float)dec[Nobj], (float)red(Nobj), (int)type[Nobj], (float)priority[Nobj],
 * (int)Nobs[Nobj].
 */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

/**< Maximum number of objects in random catalogue */
const int MAXRAND = 50000000;
/**< Maximum number of objects in mask file */
const int MAXMASK = 400000000;

int 
main (int argc, char ** argv) 
{

    FILE *maskfile;
    FILE *catalogue_file;
    FILE *random_file;
    char dd_string[1000];
    char rfilename[100];
    float *ram, *decm;
    float *rat, *dect, *zt, *prt;
    float *rar, *decr, *zr, *prr;
    int *idr, *Nobsr;
    int *idt, *Nobst;
    float dd;
    int N;
    int Nmocks = atoi(argv[4]);
    int Nmask;
    int Ncatalogue;
    int ind;
    int target_type = atoi(argv[5]);
    /**< Arrays to hold mask info */
    ram = (float *)malloc(sizeof(float)*MAXMASK);
    decm = (float *)malloc(sizeof(float)*MAXMASK);

    /**< Arrays to hold random info */
    rar = (float *)malloc(sizeof(float)*MAXRAND);
    decr = (float *)malloc(sizeof(float)*MAXRAND);
    zr = (float *)malloc(sizeof(float)*MAXRAND);
    idr = (int *)malloc(sizeof(int)*MAXRAND);
    prr = (float *)malloc(sizeof(float)*MAXRAND);
    Nobsr = (int *)malloc(sizeof(int)*MAXRAND);

     /**< Arrays to hold target info */
    rat = (float *)malloc(sizeof(float)*MAXRAND);
    dect = (float *)malloc(sizeof(float)*MAXRAND);
    zt = (float *)malloc(sizeof(float)*MAXRAND);
    idt = (int *)malloc(sizeof(int)*MAXRAND);
    prt = (float *)malloc(sizeof(float)*MAXRAND);
    Nobst = (int *)malloc(sizeof(int)*MAXRAND);
   
    maskfile = fopen(argv[1], "r");
    catalogue_file = fopen(argv[2], "r");
    printf("RAND_MAX = %d\n", RAND_MAX);

    /**< Read mask file */
    for (int i = 0; i < 1; i++) 
        fgets(dd_string, 1000, maskfile);
    N = 0;
    while (fscanf(maskfile, "%e %e %e\n", ram+N, decm+N, &dd) != EOF) 
        N++;
    Nmask = N;
    printf("read the mask file\n");

    /**< Read the target catalogue */
    N = 0;
    while (fscanf(catalogue_file, "%e %e %e %d %e %d\n", 
                  rat+N, dect+N, zt+N, idt+N, prt+N, Nobst+N) != EOF) 
        N ++;
    Ncatalogue = N;
    printf("read the catalogue file\n");

    /**< Make randoms */
    for (int i = 0; i < Nmocks; i++) {
        for (int j = 0; j < Ncatalogue; j++) {
            /**< Put all other target types in their original position */
            if (idt[j] != target_type) {
                rar[j] = rat[j];
                decr[j] = dect[j];
            }
            /** 
             * Assign random angular position from a mask to targets of this
             * type
             */
            else {
                /**< Just to make sure the selection is truly uniform */
                ind = ((double)rand())/RAND_MAX*Nmask;
                rar[j] = ram[(int)floor(ind)];
                decr[j] = decm[(int)floor(ind)];
            }
            zr[j] = zt[j];
            idr[j] = idt[j];
            prr[j] = prt[j];
            Nobsr[j] = Nobst[j];
       }
        sprintf(rfilename, "%s%02d", argv[3], i);
        random_file = fopen(rfilename, "wb");
        fwrite(&Ncatalogue, sizeof(int), 1, random_file);
        fwrite(rar, sizeof(float), Ncatalogue, random_file);
        fwrite(decr, sizeof(float), Ncatalogue, random_file);
        fwrite(zr, sizeof(float), Ncatalogue, random_file);
        fwrite(idr, sizeof(int), Ncatalogue, random_file);
        fwrite(prr, sizeof(float), Ncatalogue, random_file);
        fwrite(Nobsr, sizeof(int), Ncatalogue, random_file);
        fclose(random_file);     
        printf("created random file #%d\n", i);
    }
  
    fclose(maskfile);
    fclose(catalogue_file);

    return 0;

}
