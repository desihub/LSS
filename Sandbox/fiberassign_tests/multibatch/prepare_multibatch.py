import multibatch as mb

make_targets = False
if make_targets:
    global_DR8_mtl_file_dark = mb.make_global_DR8_mtl(output_path='targets', program='dark')
    global_DR8_mtl_file_bright = mb.make_global_DR8_mtl(output_path='targets', program='bright')
    global_DR8_sky_file = mb.make_global_DR8_sky(output_path='targets')
    global_DR8_truth_file_dark = mb.make_global_DR8_truth(global_DR8_mtl_file_dark, output_path='targets', program='dark')
    global_DR8_truth_file_bright = mb.make_global_DR8_truth(global_DR8_mtl_file_bright, output_path='targets', program='bright')

    patch_DR8_mtl_file_dark = mb.make_patch_file(global_DR8_mtl_file_dark)
    patch_DR8_mtl_file_bright = mb.make_patch_file(global_DR8_mtl_file_bright)
    patch_DR8_sky_file = mb.make_patch_file(global_DR8_sky_file)
    patch_DR8_truth_file_dark = mb.make_patch_file(global_DR8_truth_file_dark)
    patch_DR8_truth_file_bright = mb.make_patch_file(global_DR8_truth_file_bright)


make_tiles = True
if make_tiles:
    surveysim_file="/global/cfs/cdirs/desi/users/schlafly/surveysim/exposures_nopass7.fits"

    # batches for the first year of the survey (all the footprint is available) with different cadences
    n = mb.prepare_tile_batches(surveysim_file, output_path='footprint_month', program='dark', start_day=0, end_day=365, batch_cadence=28) 
    n = mb.prepare_tile_batches(surveysim_file, output_path='footprint_month', program='bright', start_day=0, end_day=365, batch_cadence=28) 
    n = mb.prepare_tile_batches(surveysim_file, output_path='footprint_semester', program='dark', start_day=0, end_day=365, batch_cadence=180) 
    n = mb.prepare_tile_batches(surveysim_file, output_path='footprint_semester', program='bright', start_day=0, end_day=365, batch_cadence=180)

    # batches for the whole duration of the survey, restricted to a small region on the sky.
    n = mb.prepare_tile_batches(surveysim_file, output_path='footprint_patch_month', program='dark', 
                                 start_day=0, end_day=2000, batch_cadence=28,select_subset_sky=True, 
                                 ra_min=130, ra_max=190, dec_min=-5, dec_max=15) 
    n = mb.prepare_tile_batches(surveysim_file, output_path='footprint_patch_month', program='bright', 
                                 start_day=0, end_day=2000, batch_cadence=28, select_subset_sky=True,
                                ra_min=130, ra_max=190, dec_min=-5, dec_max=15) 
    n= mb.prepare_tile_batches(surveysim_file, output_path='footprint_patch_semester', program='dark', 
                                 start_day=0, end_day=2000, batch_cadence=180, select_subset_sky=True,
                                ra_min=130, ra_max=190, dec_min=-5, dec_max=15) 
    n= mb.prepare_tile_batches(surveysim_file, output_path='footprint_patch_semester', program='bright', 
                                 start_day=0, end_day=2000, batch_cadence=180, select_subset_sky=True,
                               ra_min=130, ra_max=190, dec_min=-5, dec_max=15)
