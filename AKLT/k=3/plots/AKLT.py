#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')

fig, axs = plt.subplots(1, 1)

LList = [6, 8, 10, 12, 14, 16, 24, 32, 48, 64, 96, 128, 192, 256]

#------------------------------------------------------------------------------
# exact EE from MPS calculation
#------------------------------------------------------------------------------
filename = "../data/entropyMaxList.dat"
filename = "../data/entropyMaxITensor.dat"
datas = np.loadtxt(filename)
xdat, ydat, zdat = datas[:, 0], datas[:, 1], datas[:, 2]
entropyMaxList = []
for iL, L in enumerate(xdat):
    if int(L) in LList:
        entropyMaxList.append(zdat[iL])
plt.plot(LList, entropyMaxList, '-o', ms=5, label=r'$S_v^{\text{max}}$')

#------------------------------------------------------------------------------
# EE_GD bound
#------------------------------------------------------------------------------
filename = "../data/degeneracy.dat"
datas = np.loadtxt(filename)
xdat, ydat, zdat = datas[:, 0], datas[:, 1], datas[:, 2]

logD_OBCList = []
for iL, L in enumerate(xdat):
    if int(L) in LList:
        logD_OBCList.append(np.log(ydat[iL]))

plt.plot(LList, logD_OBCList, '-o', ms=5, label=r"$\log D^{\text{A}}$")

plt.xlabel(r"$L_A$")
plt.ylabel(r"$S_v$")
plt.xscale('log')
plt.legend(loc=0)
# plt.savefig('./AKLT.pdf')
plt.show()
