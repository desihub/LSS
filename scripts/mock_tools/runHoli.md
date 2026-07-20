# 1) Prepare Holi mocks

 It runs through 1000 realizations, and it will create the prepared catalogs in Y3 footprint for those available. If the catalog has been created, it won't create again.

```console
script_prepare_holi.sh
```


# 3) Run brickmask
 on all noimaging catalogs.  Run with cosmodesi enviroment branch dr1   !!!

## 3.1) Create input list
based on the catalogs we have prepared in advance, for example with: 

```console
ls /global/cfs/cdirs/desi/mocks/cai/holi/v4.00/see*/QSO/forFA0_Y3_noimagingmask_applied.fits > list_qso.txt. 
```

**Be aware many realizations have been done already for LRG, QSO. In ELG.**

## 3.2) split into batches of 70 files 
if needed

```console
split -l 70 list_qso.txt qso_list_
```

## 3.3) Create output

```console
sed 's/forFA/imforFA/g' qso_list_aa > out_list_aa.txt
```

## 3.4) Execute brickmask


```console
brickmask/script_holi.sbatch

```
 with

```console
config brickmask_holi.conf
```

maybe several times depending on the number of lists. Need to change INPUT_FILES and OUTPUT_FILES command to the list of files

**Remeber to change parameters INPUT_FILES and OUTPUT_FILES!!!**




# 4) Apply NOBS, MASKBIT mask

> NOTE
>
> run 4,5,6 with combine_step456_z1.sl


Apply NOBS, MASKBIT mask on files created in step 2

```console
python join_imaging_mask.py
```

You need to change the path to where to find the list of files, defined in OUTPUT_FILES in the previous step. You can put several files in a list, and it will run through all of them. 

For example, in step 3, you needed to define an output list for ELG, another for LRG and another for QSO (maybe even more!). This will run on every file in the lists, and save it with name forFA


# 5) Add contaminants 

to ELG and QSO samples. It runs over the realizations shown in the discovery step

```console
python add_contaminants_to_mock.py
```


# 6) Concatenate tracers
 into a forFA file (saving to scratch and creating the qso file needed in AltMTL)

```console
python concatenate_tracers_to_fba.py
```

> NOTE
>
>run 4,5,6 with combine_step456_z1.sl

# 7) Initialize the altmtl directories 

this will done over the whole list discovered previously. Or if you want, just select the first 100 that w

```console
nohup ./script_init_holi.sh &
```

It's a loop on script for each seed catalog

```console
initialize_amtl_mocks_da2.py
```

# 8) Run altMTL

## Recent script

```console
scripts/mock_tools/DR2_altmtl_sbatch/altmtl_200249_holi_v3.sbatch
``` 
[github link](https://github.com/desihub/LSS/blob/main/scripts/mock_tools/DR2_altmtl_sbatch/altmtl_200249_holi_v3.sbatch)

## Previous script

```console
DA2ALTMTLRealizationsDARK_holi.sh
``` 
 
 in LSS/bin  The list of 100 realizations is already preprated. Just execute it normaly. It will save to scratch directory

## chatGPT Analysis of LSS/bin/DA2ALTMTLRealizationsDARK_holi.sh.

Main 4 processing stages
1. Environment and machine resource setup  
   The script loads the DESI environment, detects the NERSC machine (cori/perlmutter), sets batch parameters (constraints, queue, ProcPerNode), then initializes debug/verbose/mock flags.

2. Realization campaign parameterization  
   It defines the mock range (mockinit/mockend or mocklist), observing condition (DARK/BRIGHT), survey (main), shuffle/subpriority options, restart options, and automatically computes the number of nodes NNodes.

3. Output directory construction and reproducibility handling  
   It builds output paths from ALTMTLHOME + simName, creates directories, copies the script into the output folder, and checks consistency with an already copied version (cmp) before continuing.

4. AltMTL production launch  
   It assembles argstring (including targfile, zfix, mock range, etc.) and launches the pipeline through dateLoopAltMTLBugFix_mock_batch.sh using nohup, then performs status/runtime checks.

Paths to change to run on your own simulated data
In LSS/bin/DA2ALTMTLRealizationsDARK_holi.sh, the most important ones are:

1. LSS code location  
   path2LSS=/global/homes/d/desica/LSScode/LSS/bin/  
   Replace with your local clone path.

2. AltMTL output root directory  
   ALTMTLHOME=/pscratch/sd/d/desica/DA2/mocks/holi_v3/  
   Replace with your production output directory.

3. Mock target input forFA  
   targfile="--targfile=/pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA{mock_number}.fits"  
   Replace with the path to your forFAxxxx.fits files.

4. QSO z-fix input (optional but often required)  
   zfix="--zfix=/pscratch/sd/d/desica/DA2/mocks/holi_v3/qsos/qso{mock_number}.txt"  
   Replace with your qsos directory.

5. Campaign/logical naming  
   simName="altmtl{mock_number}"  
   Adapt if you want a different naming convention.

6. DESI environment source (if different on your side)  
   source desi_environment.sh main  
   Keep if valid for you, otherwise adapt.

7. Mock range  
   mockinit=0, mockend=50 or mocklist  
   Adjust to your available realizations.

Important caveats
1. The placeholder {mock_number} in simName, targfile, zfix is not replaced by this script as shown.  
   You need either pre-substitution before execution, or an explicit shell variable and string construction with that variable.

2. This line is likely incorrect:  
   cp $SLURM_SUBMIT_DIR $0 $outputMTLFinalDestination/$0  
   It passes 2 sources to cp. Usually you only want to copy the script itself.

3. endInit is used in runtime calculations but is not defined in the shown script.  
   That can break timing logic.

4. The script ends with exit 54321 even on success, which is a non-zero exit code.  
   In batch systems, this may be interpreted as failure.


# 9) In parallel, create pota and lrg mask files

put realizations accordingly: 

```console
 /global/homes/d/desica/LSScode/LSS/scripts/mock_tools/lrgmask_pota_holi.sh
```

using scripts:

```console
abamtl_getpota_sbatch_da2_holi.sh abamtl_lrgmask_sbatch_da2_holi.sh
```

AFTER IT FINISHES

# 10) Eliminates spurious temporary directories

remove 

```console
altmtlX/initled
holi/tartilesX
```


NEXT: RUN LSS PIPELINE (DIFFERENT SCRIPT)
