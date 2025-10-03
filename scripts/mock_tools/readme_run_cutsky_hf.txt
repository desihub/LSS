python prep_cutsky_torun.py --input_mockpath /global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0 --input_mockfile forFA0_Y3_noimagingmask_applied.fits --tracer LRG --snapshot z0.500 --abacus_realization 0 --redshift_range 0.4,1.1 --nrans 4




python prepare_mocks_Y3.py --mockname abacushf --input_mockpath /global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO/z1.400 --input_mockfile cutsky_QSO_z1.400_AbacusSummit_base_c000_ph000.fits --output_fullpathfn /global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO/z1.400/newtest.fits --tracer QSO --zrsdcol Z --ztruecol Z_COSMO --need_footprint n --need_nz_calib n


srun -N 1 -C cpu --cpus-per-task=128 -t 04:00:00 --ntasks=1 --qos interactive --account desi python apply_phot_mask.py 128 1 DARK


python concatenate_tracers_to_fba.py


#modified by hand, hardwire 
srun -N 1 -C cpu --cpus-per-task=256 -t 04:00:00 --ntasks=1 --qos interactive --account desi python readwrite_pixel_bitmask_da2.py --tracer lrg -i 0 --cat_type Ab2ndgen --secgen_ver AbacusSummit_v4_1 --nproc 256


srun -N 1 -C cpu --cpus-per-task=64 -t 04:00:00 --ntasks=1 --qos interactive --account desi python initialize_amtl_mocks_da2.py /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0_withQSOcont.fits /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/altmtl0_withQSOcont DARK
