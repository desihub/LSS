import sys
import os
sys.path.append(os.path.abspath('multibatch'))
import multibatch as mb

target_path='targets'
surveysim_file="/global/cfs/cdirs/desi/users/schlafly/surveysim/exposures_nopass7.fits"
bdir = ''
ra_min=130
ra_max=190
dec_min=-5
dec_max=15
batch_cadence=28
if batch_cadence == 28:
   foot =bdir+'footprint_month'
   darkout = bdir+"test_dark_global_month"

make_glob = False
if make_glob:
    global_DR8_mtl_file_dark = mb.make_global_DR8_mtl(output_path=target_path, program='dark')
    global_DR8_mtl_file_bright = mb.make_global_DR8_mtl(output_path=target_path, program='bright')    
    global_DR8_sky_file = mb.make_global_DR8_sky(output_path=target_path)
    global_DR8_truth_file_dark = mb.make_global_DR8_truth(global_DR8_mtl_file_dark, output_path=target_path, program='dark')
    global_DR8_truth_file_bright = mb.make_global_DR8_truth(global_DR8_mtl_file_bright, output_path=target_path, program='bright')
else:
    output_path=target_path 
    program='dark'
    global_DR8_mtl_file_dark = os.path.join(output_path, 'global_DR8_mtl_{}.fits'.format(program))
    global_DR8_truth_file_dark = os.path.join(output_path, 'global_DR8_truth_{}.fits'.format(program))
    program = 'bright'
    global_DR8_mtl_file_bright = os.path.join(output_path, 'global_DR8_mtl_{}.fits'.format(program))
    global_DR8_truth_file_bright = os.path.join(output_path, 'global_DR8_truth_{}.fits'.format(program))
    global_DR8_sky_file = os.path.join(output_path, "global_DR8_sky.fits")
     
make_targets = False
if make_targets:

    patch_DR8_mtl_file_dark = mb.make_patch_file(global_DR8_mtl_file_dark,ra_min, ra_max, dec_min, dec_max)
    patch_DR8_mtl_file_bright = mb.make_patch_file(global_DR8_mtl_file_bright,ra_min, ra_max, dec_min, dec_max)
    patch_DR8_sky_file = mb.make_patch_file(global_DR8_sky_file,ra_min, ra_max, dec_min, dec_max)
    patch_DR8_truth_file_dark = mb.make_patch_file(global_DR8_truth_file_dark,ra_min, ra_max, dec_min, dec_max)
    patch_DR8_truth_file_bright = mb.make_patch_file(global_DR8_truth_file_bright,ra_min, ra_max, dec_min, dec_max)


make_tiles = False
if make_tiles:
    

    # batches for the first year of the survey (all the footprint is available) with different cadences
    n = mb.prepare_tile_batches(surveysim_file, output_path=foot, program='dark', start_day=365, end_day=730, batch_cadence=batch_cadence) 
    #n = mb.prepare_tile_batches(surveysim_file, output_path=foot, program='bright', start_day=365, end_day=730, batch_cadence=batch_cadence) 

    # batches for the whole duration of the survey, restricted to a small region on the sky.
    #n = mb.prepare_tile_batches(surveysim_file, output_path='footprint_patch_month', program='dark', 
    #                             start_day=0, end_day=2000, batch_cadence=batch_cadence,select_subset_sky=True, 
    #                            ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max) 
    #n = mb.prepare_tile_batches(surveysim_file, output_path='footprint_patch_month', program='bright', 
    #                             start_day=0, end_day=2000, batch_cadence=batch_cadence, select_subset_sky=True,
    #                            ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max) 

run_fa = True
if run_fa:
    sbatch = int(sys.argv[1])
    mxbatch = int(sys.argv[2])
    #for i in range(sbatch,mxbatch+1):
    mb.run_strategy(global_DR8_mtl_file_dark,  global_DR8_truth_file_dark , global_DR8_sky_file,
        output_path=darkout, batch_path=foot, program="dark",sbatch=sbatch,mxbatch=mxbatch)#,sbatch=i,mxbatch=i+1)
