#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This program plots Fig. S2 in the supplementary material."""
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')
plt.rcParams.update({
    "font.family": "sans-serif",     # specify font family here
    "font.sans-serif": ["Arial"]})   # specify font here
plt.rcParams['figure.autolayout'] = False

#--------------------------------------------------------------------------
L, k = 20, 5
iscarList = [1, 2, 7, 30, 102, 354, 426, 449, 454, 455]
#--------------------------------------------------------------------------

filename   = f"../data/L{L}_k{k}_lamMat.dat"
lamvalMat  = np.loadtxt(filename)
dimLH, basisN = np.shape(lamvalMat)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.8, 5.2))

def plot(ax, n, ymax=0):
    """
    plot correlation matrix spectrum for few eigenlevels around n
    """
    iList = []
    for j in range(-3, 3+1, 1):
        i = n + j
        if i in range(basisN):
            iList.append(i)
            if i in iscarList:
                ax.plot(list(range(1, len(lamvalMat[:, i-1])+1)), lamvalMat[:, i-1], 'ko-', ms=2, lw=1, label=r"$n="+str(i)+"$")
            else:
                ax.plot(list(range(1, len(lamvalMat[:, i-1])+1)), lamvalMat[:, i-1], 'o-', ms=2, label=r"$n="+str(i)+"$")
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\lambda_\alpha$')
    ax.legend(loc=0, ncol=2)
    if ymax == 0:
        ymax = max(lamvalMat[-1, iList])
    ax.set_ylim(top=ymax, bottom=-0.5)

plot(ax1, 30, 6)
plot(ax2, 102, 6)

ax1.text(-0.12, 1.02, '(a)', transform=ax1.transAxes, fontsize = '12')
ax1.text(-0.12, 1.02, '(b)', transform=ax2.transAxes, fontsize = '12')

fig.subplots_adjust(hspace=0.25)
# plt.savefig(f"./PXP_L{L}_k{k}.pdf")
plt.show()
