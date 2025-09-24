import multibatch as mb


mb.run_strategy("targets/global_DR8_mtl_dark.fits", "targets/global_DR8_truth_dark.fits", "targets/global_DR8_sky.fits",
             output_path="test_dark_global_semester", batch_path="footprint_semester", program="dark")

#mb.run_strategy("targets/global_DR8_mtl_dark.fits", "targets/global_DR8_truth_dark.fits", "targets/global_DR8_sky.fits", 
#             output_path="test_dark_global_month", batch_path="footprint_month", program="dark")

#mb.run_strategy("targets/patch_DR8_mtl_dark.fits", "targets/patch_DR8_truth_dark.fits", "targets/patch_DR8_sky.fits",
#             output_path="test_dark_patch_semester", batch_path="footprint_patch_semester", program="dark")

#mb.run_strategy("targets/patch_DR8_mtl_dark.fits", "targets/patch_DR8_truth_dark.fits", "targets/patch_DR8_sky.fits", 
#             output_path="test_dark_patch_month", batch_path="footprint_patch_month", program="dark")

