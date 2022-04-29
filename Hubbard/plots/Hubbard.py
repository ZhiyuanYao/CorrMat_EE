#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')

fig, axs = plt.subplots(1, 1)

LList = [6, 8, 10, 12, 14, 16, 24, 32, 48, 64, 96, 128, 192, 256]

#------------------------------------------------------------------------------
# ground-state degeneracy data
#------------------------------------------------------------------------------
datas = np.loadtxt("../data/Hubbard.dat")
xdat, ydat, zdat = datas[:, 0], datas[:, 1], datas[:, 2]
entropyList = []
logDOBCList = []
for iL, L in enumerate(xdat):
    if int(L) in LList:
        entropyList.append(ydat[iL])
        logDOBCList.append(zdat[iL])

plt.plot(LList, entropyList, '-o', ms=5, label=r'$S_v$')
plt.plot(LList, logDOBCList, '-o', ms=5, label=r'$\log D^{\text{A}}$')

plt.xlabel(r"$L_A$")
plt.ylabel(r"$S_v$")
plt.xscale('log')
plt.legend(loc=0)
# plt.savefig('Hubbard.pdf')
plt.show()
