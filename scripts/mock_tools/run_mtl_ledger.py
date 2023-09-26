from desitarget import mtl
import sys
import glob

arg1 = sys.argv[1] #Input mock
arg2 = sys.argv[2] #Output path
obscon = sys.argv[3] #DARK or BRIGHT

print('Running initial ledgers')

mtl.make_ledger(arg1, arg2, obscon=obscon.upper())

print('Creating list of tiles to be processed by AltMTL mock production')

path = os.path.join(arg2, 'main', obscon.lower())

ff = glob.glob(os.path.join(path, 'mtl-{obscon}-hp-*.ecsv'.format(obscon=obscon.lower())))

dd=[]

for f in ff:
    dd.append(int(f.split('hp-')[-1].split('.ecsv')[0]))
tosave = ','.join(map(str, sorted(dd)))

savepath = os.path.join(arg2, 'hpxlist.txt')

ff = open(savepath, 'w')
ff.write(tosave)
ff.close()

print('saving list of HP ledgers in '+savepath)
