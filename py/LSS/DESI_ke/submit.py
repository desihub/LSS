import os
import argparse
import numpy as np


def customise_script(args, debug=False):
    script     = args.script
    queue      = args.queue
    memory     = args.memory
    time       = args.time
    script_log = args.script_log
    account    = args.account
    nodes      = args.nodes

    root       = os.environ['CODE_ROOT'] + '/bin/'

    supported  = ['gold_pipeline',\
                  'gold_d8_pipeline',\
                  'rand_pipeline',\
                  'rand_d8_pipeline',\
                  'rand_ddp1_d8_pipeline',\
                  'rand_ddp1_pipeline']

    if script == None:
        for script in supported:
            args.script = root + script

            customise_script(args)

        return 0

    else:
        assert root in script, f'Error on {script}'

    if script_log == None:
        #  TODO:  findfile                                                                                                                                                                               
        ss   = os.path.basename(script).split('.')[0]

        if 'rand' in script:
            log  = os.environ['HOME'] + f'/data/GAMA4/randoms/logs/{ss}.log'
        else:
            log  = os.environ['HOME'] + f'/data/GAMA4/logs/{ss}.log'
    else:
        log  = script_log

    ## 
    ff       = open(f'/{script}')
    ff       = ff.read()
    ff       = ff.split('\n')
    ff       = [x.rstrip() for x in ff]

    custom   = [x for x in ff if '#SBATCH' in x]
    rest     = [x for x in ff if '#SBATCH' not in x] 
    
    #
    #  --  Example  --
    #
    #SBATCH -p cordelia
    #SBATCH --mem=20G
    #SBATCH -t 02:00:00
    #SBATCH -o /cosma/home/durham/dc-wils7/data/GAMA4/logs/gold_pipeline.log
    #SBATCH -A durham
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --open-mode=append
    
    def swap(xx):
        args = {'queue':   {'sbatch': ['-p', '--partition'], 'arg': queue,        'equals': False},
                'memory':  {'sbatch': ['--mem'],             'arg': memory,       'equals': True},
                'time':    {'sbatch': ['-t'],                'arg': time,         'equals': False},
                'account': {'sbatch': ['-A'],                'arg': account,      'equals': False},
                'log':     {'sbatch': ['-o'],                'arg': log,          'equals': False},
                'nodes':   {'sbatch': ['--nodes'],           'arg': nodes,        'equals': True}}
        
        if debug:
            print('\n\n...  Swapping  ...')

        for i, _ in enumerate(xx):
            for arg in sorted(args.keys()):
                var    = args[arg]['arg']
                batch  = np.atleast_1d(args[arg]['sbatch'])
                equals = args[arg]['equals'] 

                if (var != None):
                    if arg == 'memory':
                        var += 'G'

                    split = []

                    for yy in _.split():
                        split += yy.split('=')

                    found = np.any([xx in batch for xx in split])
                    
                    # print(var, batch, found, split)

                    if found:
                        if equals:
                            xx[i] = '#SBATCH {}={}'.format(batch[0], var)

                        else:
                            xx[i] = '#SBATCH {} {}'.format(batch[0], var)

        return  xx

    if debug:
        print('\n\n----  INPUT  ----\n')

        for xx in custom:
            print(xx)


    ##  Swap
    custom = swap(custom)

    if debug:
        print('\n\n----  OUTPUT  ----\n')

        for xx in custom:
            print(xx)


    opath = root + '/custom/' + os.path.basename(script).split('.')[0]
            
    with open(opath, 'w') as f:
        rest.remove('#!/bin/bash')
        
        to_write = custom + rest

        f.write('#!/bin/bash')
        f.write('\n')

        for line in custom + rest:
            f.write(line)
            f.write('\n')

    print('\nWriting {}'.format(opath))

    os.system(f'chmod 700 {opath}')


if __name__ == '__main__':
    # /cosma/home/durham/dc-wils7/DESI/bin/gold_pipeline                                                                                                                                                 
    parser    = argparse.ArgumentParser(description='Customise pipeline submission scripts.')
    parser.add_argument('-s', '--script',  help='Script to customise.',    type=str, default=None)
    parser.add_argument('-q', '--queue',   help='Queue for submission.',   type=str, default=None)
    parser.add_argument('-m', '--memory',  help='Node memory usage [GB].', type=str, default=None)
    parser.add_argument('-t', '--time',    help='Job time to request.',    type=str, default=None)
    parser.add_argument('--script_log',    help='Job log path.',           type=str, default=None)
    parser.add_argument('-a', '--account', help='Account for submission.', type=str, default=None)
    parser.add_argument('-n', '--nodes',   help='Nodes to request.',       type=int, default=None)

    args      = parser.parse_args()
    script    = args.script

    customise_script(args)

    print('\n\nDone.\n\n')
