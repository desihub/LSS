import numpy as np

from   delta8_limits import d8_limits


def volfracs(rand, bitmasks=['IN_D8LUMFN']):
    utiers    = np.unique(rand['DDP1_DELTA8_TIER'].data)
    utiers_zp = np.unique(rand['DDP1_DELTA8_TIER_ZEROPOINT'].data)

    utiers    = utiers[utiers >= 0]

    print('Unique tiers: {}'.format(utiers))

    ddp1_rand = rand[rand['DDPZLIMS'][:,0] == 1]

    print('DDP1 randoms: {:.6f} < z < {:.6f}'.format(ddp1_rand['Z'].min(), ddp1_rand['Z'].max()))
    
    for ut in range(len(d8_limits)):
        print()

        # TODO/HACK/BUG/MJW  Shouldn't this be ddp1_rand?  Rather than rand.
        in_tier = (ddp1_rand['DDP1_DELTA8_TIER'].data == ut)

        for bm in bitmasks:
            in_tier &= (ddp1_rand[bm].data == 0)

            print(bm, np.mean(in_tier))

        # print(ut, d8_limits[ut], np.mean(d8_limits[ut]))

        if np.count_nonzero(in_tier) > 0: 
            rand.meta['DDP1_d{}_VOLFRAC'.format(ut)]   = '{:.6f}'.format(np.mean(in_tier))
            rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)] = '{:.6f}'.format(np.mean(ddp1_rand['DDP1_DELTA8'].data[in_tier]))

        else:
            rand.meta['DDP1_d{}_VOLFRAC'.format(ut)]   = '{:.6f}'.format(0.0)
            rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)] = '{:.6f}'.format(np.mean(d8_limits[ut]))
            
        print('DDP1_d{}_VOLFRAC OF {} added.'.format(ut,    rand.meta['DDP1_d{}_VOLFRAC'.format(ut)]))
        print('DDP1_d{}_TIERMED d8 OF {} added.'.format(ut, rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)]))

        # Zero point.                                                                                                                                                                                      
        in_tier = (ddp1_rand['DDP1_DELTA8_TIER_ZEROPOINT'].data == ut) 

        for bm in bitmasks:
            in_tier &= (ddp1_rand[bm].data == 0)

        if np.count_nonzero(in_tier) > 0:
            rand.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(ut)]   = '{:.6f}'.format(np.mean(in_tier))
            rand.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(ut)] = '{:.6}'.format(np.mean(ddp1_rand['DDP1_DELTA8_ZEROPOINT'].data[in_tier]))
        
        else:
            rand.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(ut)]   = '{:.6f}'.format(0.0)
            rand.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(ut)] = '{:.6f}'.format(np.mean(d8_limits[ut]))

        print('DDP1_d{}_ZEROPOINT_VOLFRAC OF {:.10f} added.'.format(ut, np.mean(in_tier)))
        print('DDP1_d{}_ZEROPOINT_TIERMED d8 OF {} added.'.format(ut, rand.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(ut)]))

    return  rand
