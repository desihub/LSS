# A code to compute correlation function moments xi0 and xi2 from pair-count
# files for DD=argv[1], DR=argv[2] and RR=argv[3]. The output xi0 file = argv[4]
# and the output xi2 file = argv[5]

from sys import argv
import numpy as np
import string as str

Nbins = 30

DDfile = argv[1]
DRfile = argv[2]
RRfile = argv[3]
xi0file = argv[4]
xi2file = argv[5]

DDf = open(DDfile)
header = DDf.readline()
header = header.split()
NDD = str.atof(header[3])
header = DDf.readline()
header = header.split()
NR = str.atof(header[2])
header = DDf.readline()
header = header.split()
NM = str.atof(header[2])
DDf.close()

DRf = open(DRfile)
header = DRf.readline()
header = header.split()
NDR = str.atof(header[3])
header = DRf.readline()
header = DRf.readline()
DRf.close()

RRf = open(RRfile)
header = RRf.readline()
header = header.split()
NRR = str.atof(header[3])
header = RRf.readline()
header = RRf.readline()
RRf.close()

DD = np.loadtxt(argv[1])
DR = np.loadtxt(argv[2])
RR = np.loadtxt(argv[3])

print NDD, NDR, RR.shape
DDnew = np.zeros((NM, Nbins))
DRnew = np.zeros((NM, Nbins))
RRnew = np.zeros((NM, Nbins))

for i in range(0, Nbins):
    for j in range(0, int(NM)):
        DDnew[j,i] = sum(DD[j,i*5:(i+1)*5])
        DRnew[j,i] = sum(DR[j,i*5:(i+1)*5])
        RRnew[j,i] = sum(RR[j,i*5:(i+1)*5])

xi2d = (NRR/NDD*DDnew - 2*NRR/NDR*DRnew + RRnew)/RRnew
print xi2d

mu = np.linspace(1.0/NM/2.0, 1.0 - 1.0/NM/2.0, NM)
L2 = (3.0*mu*mu - 1.0)/2.0
L0 = mu + 1.0 - mu;

xi0 = np.dot(L0, xi2d)/NM
print xi0.shape
print xi0

np.savetxt(xi0file, xi0)

xi2 = np.dot(L2, xi2d)/NM*5
print xi2.shape
print xi2

np.savetxt(xi2file, xi2)

