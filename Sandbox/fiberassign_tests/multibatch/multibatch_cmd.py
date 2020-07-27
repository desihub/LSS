from __future__ import absolute_import, division, print_function


import os,argparse,time
import multibatch as mb
parser = argparse.ArgumentParser()

# Required arguments
parser.add_argument("--mtl",help="input targets (FITS file)",type=str, required=True)
parser.add_argument("--sky",help="input sky positions (FITS file)",type=str, required=True)
parser.add_argument("--truth",help="input targets truth (FITS file)",type=str, required=True)
parser.add_argument("--expfile",type=str,help="input exposures file path (FITS file)", required=True)

# Optional arguments
parser.add_argument("--outdir",help="output directory",type=str,default='./')
parser.add_argument("--cadence",help="Cadence for multiple batches generation and fiber assignment",type=int, default=None)
parser.add_argument("--first-day",help="Starting day of the period for fiber assignment",type=int, default=0)
parser.add_argument("--last-day",help="Last day of the period for fiber assignment",type=int, default=365)
parser.add_argument("--program",help="Program name, e.g dark or bright",type=str, default='dark')

args = parser.parse_args()

print('Starting at {}'.format(time.asctime()))
exp_outdir=os.path.join(args.outdir,'exposures')
print("Saving batches files in {}".format(exp_outdir))
n = mb.prepare_tile_batches(surveysim_file=args.expfile, output_path=exp_outdir, 
                            start_day=args.first_day, end_day=args.last_day, batch_cadence=args.cadence)
print("Running multibatch")
mb.run_strategy(args.mtl, args.truth, args.sky, output_path=args.outdir, batch_path=exp_outdir, program=args.program.lower())
print('Finished at {}'.format(time.asctime()))