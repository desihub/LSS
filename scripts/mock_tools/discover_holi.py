import os

all_realizations = list(range(1000))


parent_path = '/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80'
#holi_path = '/global/cfs/cdirs/desi/mocks/cai/holi/'


vQSO = 'v4.00'
vELG = 'v4.00'
vLRG = 'v4.00'

todo = []

for i in all_realizations:
    seednum = str(i).zfill(4)
    tempdir = os.path.join(parent_path, 'seed%s' % seednum)
    #print(tempdir)
    if os.path.isdir(tempdir):
        if os.path.isfile(os.path.join(tempdir, 'holi_LRG_v4.80_GCcomb_clustering.dat.h5')) and os.path.isfile(os.path.join(tempdir, 'holi_QSO_v4.80_GCcomb_clustering.dat.h5')) and os.path.isfile(os.path.join(tempdir, 'holi_ELG_v4.80_GCcomb_clustering.dat.h5')):
            #todo.append(i)
    #if not os.path.isdir(os.path.join(parent_path, 'altmtl%d' % i)):
#        if os.path.isdir(os.path.join(holi_path, vQSO, 'seed%s' % seednum)) and os.path.isdir(os.path.join(holi_path, vELG, 'seed%s' % seednum)):
            todo.append(i)

print(todo)
print(len(todo))

strr = ''
for i in todo:
    strr += '%d '%i

print(strr)

