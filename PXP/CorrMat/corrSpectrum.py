#!/usr/bin/env python
"""This program plots the corr matx spectrum with black line for scars."""
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')
# plt.rc('text', usetex=True)

#--------------------------------------------------------------------------
L, k = 20, 5
iscarList = [1, 2, 7, 30, 102, 354, 426, 449, 454, 455]
#--------------------------------------------------------------------------

filename   = f"./data/L{L}_k{k}_lamMat.dat"
lamvalMat  = np.loadtxt(filename)
dimLH, basisN = np.shape(lamvalMat)

def plot(n, ymax=0):
    """
    plot correlation matrix spectrum for few eigenlevels around n
    """
    fig, ax = plt.subplots(1, 1, figsize=(4.8, 3.2))
    iList = []
    for j in range(-3, 3+1, 1):
        i = n + j
        if i in range(basisN):
            iList.append(i)
            if i in iscarList:
                plt.plot(list(range(1, len(lamvalMat[:, i-1])+1)), lamvalMat[:, i-1], 'ko-', ms=2, lw=1, label="n="+str(i))
            else:
                plt.plot(list(range(1, len(lamvalMat[:, i-1])+1)), lamvalMat[:, i-1], 'o-', ms=2, label=f"n={i}")
    plt.xlabel(r'$i$')
    plt.ylabel(r'$\lambda_i$')
    # plt.title(f"L = {L}, range = {k}", fontsize=10)
    plt.legend(loc=0)
    if ymax == 0:
        ymax = max(lamvalMat[-1, iList])
    plt.gca().set_ylim(top=ymax, bottom=-0.5)
    plt.savefig(f"./plots/L{L}_k{k}_{n}.pdf")
    plt.show()


def plots(nlist, ymax=0):
    """
    plot correlation matrix spectrum for eigenlevels in nlist
    """
    fig, ax = plt.subplots(1, 1, figsize=(4.8, 3.2))
    iList = []
    for i in nlist:
        iList.append(i)
        if i in iscarList:
            plt.plot(list(range(1, len(lamvalMat[:, i-1])+1)), lamvalMat[:, i-1], 'ko-', ms=2, lw=1, label="n="+str(i))
        else:
            plt.plot(list(range(1, len(lamvalMat[:, i-1])+1)), lamvalMat[:, i-1], 'o-', ms=2, label=f"n={i}")
    plt.xlabel(r'$i$')
    plt.ylabel(r'$\lambda_{i}$')
    # plt.title(f"L = {L}, range = {k}", fontsize=10)
    plt.legend(loc=0)
    if ymax == 0:
        ymax = max(lamvalMat[-1, iList])
    plt.gca().set_ylim(top=ymax, bottom=-0.5)
    plt.savefig(f"./plots/L{L}_k{k}_{nlist[0]}.pdf")
    plt.show()

# plt.savefig(./plots/L20_k5_30.pdf")
