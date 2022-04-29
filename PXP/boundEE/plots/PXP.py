#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This program is to fit and draw the (xname, yname) line."""
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=11)
plt.rc('ytick', labelsize=11)
plt.rcParams["font.family"] = "Times New Roman"

L = 20
iscarList = [1, 2, 7, 30, 102, 354, 426, 449, 454, 455]    # kI= 0+
#------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 1)
#------------------------------------------------------------------------------
filename = f"./Es_EE_L{L}_m0.000_k=0+.dat"
datas = np.loadtxt(filename)
EList, SList = datas[:, 0], datas[:, 1]
plt.plot(EList, SList, 'o', ms=4, alpha=0.6)

for k, marker in zip([3, 4, 5], ['p-', 's-', 'D-']):
    datas = np.loadtxt(f"../data/L{L}_k{k}_EE_bound.dat")
    nlist, entropyList = datas[:, 0], datas[:, 1]
    nlist = [i-1 for i in iscarList]
    entropyList = np.append(entropyList, entropyList[::-1])
    plt.plot(EList[nlist], entropyList, marker, ms=4, label=f"range={k}")
    if k == 5:
        plt.plot(EList[nlist], SList[nlist], 'kx', ms=4, label="scar states")

plt.legend(loc=2, ncol=2)
plt.xlabel(r'$E_n$')
plt.ylabel(r'$S_v$')
# plt.title("xname-yname plot")

plt.savefig(f"PXP.pdf")
# plt.show()
