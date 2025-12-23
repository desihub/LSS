import os

all_realizations = list(range(1000))


parent_path = '/global/cfs/cdirs/desi/mocks/cai/LSS/DA2/mocks/holi_v1'
holi_path = '/global/cfs/cdirs/desi/mocks/cai/holi/'


vQSO = 'v4.00'
vELG = 'v5.0'
vLRG = 'v4.00'

todo = []

for i in all_realizations:
    seednum = str(i).zfill(4)
    if not os.path.isdir(os.path.join(parent_path, 'altmtl%d' % i)):
        if os.path.isdir(os.path.join(holi_path, vQSO, 'seed%s' % seednum)) and os.path.isdir(os.path.join(holi_path, vELG, 'seed%s' % seednum)):
            todo.append(i)

print(todo)
print(len(todo))

strr = ''
for i in todo:
    strr += '%d '%i

print(strr)

