# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 10:24:56 2022

@author: Santi

Version Chevallier de la FISTA 1D para python.

Comentarios de la FISTA1D original:

% Version1D a partir de
% Fast 2D NMR relaxation distribution estimation - Matlab/octave version
% Paul Teal, Victoria University of Wellington
% paul.teal@vuw.ac.nz
% Let me know of feature requests, and if you find this algorithm does
% not perform as it should, please send me the data-set, so I can improve it.
% Issued under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3.
% If you distribute this code or a derivative work you are required to make the
% source available under the same terms
% If you use this software, please cite P.D. Teal and C. Eccles. Adaptive
% truncation of matrix decompositions and efficient estimation of NMR
% relaxation distributions. Inverse Problems, 31(4):045010, April
% 2015. http://dx.doi.org/10.1088/0266-5611/31/4/045010 (Section 4: although
% the Lipshitz constant there does not have alpha added as it should have)

% Y is the NMR data for inversion
% alpha is the (Tikhonov) regularisation (scalar)
% S is an optional starting estimate

% K1 is the kernel matrices
% They can be created with something like this:

"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler
from scipy.signal import find_peaks
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.optimize import curve_fit



# K = np.zeros((tau.size, Npts_L))  # Npts_L: puntos de la transformada
# for ii in tau.size:
#     K[ii,:] = 1 - np.exp(-tau[ii] / T)

def NLI_FISTA_1D(K, Z, alpha, S):
    '''
    Inversión de Laplace 1D
    '''

    Z = np.reshape(Z, (len(Z), 1))
    S = np.reshape(S, (len(S), 1))

    KTK = K.T @ K # @: matrix multiplication: @ = np.dot()
    KTZ = K.T @ Z
    ZZT = np.trace(Z @ Z.T)

    invL = 1 / (np.trace(KTK) + alpha)
    factor = 1 - alpha * invL

    Y = S
    tstep = 1
    lastRes = np.inf

    for iter in range(100000):
        term2 = KTZ - KTK @ Y
        Snew = factor * Y + invL * term2
        Snew[Snew<0] = 0

        tnew = 0.5 * (1 + np.sqrt(1 + 4 * tstep**2))
        tRatio = (tstep - 1) / tnew
        Y = Snew + tRatio * (Snew - S)
        tstep = tnew
        S = Snew

        if iter % 1000 == 0:
            TikhTerm = alpha * np.linalg.norm(S)**2
            ObjFunc = ZZT - 2 * np.trace(S.T @ KTZ) + np.trace(S.T @ KTK @ S) + TikhTerm

            Res = np.abs(ObjFunc - lastRes) / ObjFunc
            lastRes = ObjFunc
            print(f'# It = {iter} >>> Residue = {Res:.6f}')

            if Res < 1E-5:
                break

    return S[:, 0]

def fitMag_1D(tau1, T1, S_1D):
    '''
    Ajuste del decaimiento a partir de la distribución de T1 de la Laplace 1D.
    '''
    M = []
    for i in range(tau1.size):
        m = 0
        for j in range(T1.size):
            m += S_1D[j] * (1 - np.exp(- tau1[i] / T1[j]))
        M.append(m)

    return np.array(M)


#%%
# de antemano tengo las variables tau, signal, tau_fit


Npts_L = 512
alpha = 1e-4

T = np.logspace(1,5, Npts_L)

K = np.zeros((tau.size, Npts_L))  # Npts_L: puntos de la transformada
for ii in range(tau.size):
    K[ii,:] = (1 - np.exp(-tau[ii] / T))


S = T*0
S = NLI_FISTA_1D(K, signal, alpha, S)

s_ilt = fitMag_1D(tau_fit, T, S)

fig, axs = plt.subplots(1, 2)

axs[0].semilogx(T,S, 'r-')

axs[1].plot(tau, signal, 'o')
axs[1].plot(tau_fit, s_ilt, 'r-')


