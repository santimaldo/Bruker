import nmrglue as ng
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from Datos import *
from Exsy import *
from Espectro import autophase
from scipy import integrate
import matplotlib.cm




from scipy.signal import savgol_filter


def Integrar(Matriz, x=None, y=None):
    """
    Integracion 2d
    """
    if x is None:
        i = (integrate.simpson(integrate.simpson(Matriz)))
    else:
        i = (integrate.simpson(integrate.simpson(Matriz, x=x), x=y))
    return i

# #############################################################################
# #############################################################################
# ###################CMK3ACTB 1H  # muestra "B". con cantidad de bulk correcta
# #############################################################################
# #############################################################################
# -1 V
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-03-10_insitu-LiTFSIaq-supercap/"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\exsy/-1V/"
mT = [0.1, 1, 10, 100, 5, 50, 200, 20, 2, 80, 150]  # 1H
nexp = np.arange(80,91)
peaks = [0, -5.9]
# 0 V
# path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-03-10_insitu-LiTFSIaq-supercap/"
# savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\exsy/0V/"
# mT = [0.1, 1, 10, 100, 5, 50, 200, 20, 2, 80, 150]  # 1H
# nexp = np.arange(20,31)
# peaks = [0, -5]
semiancho = 1


# reordeno - -------------------
# zipped_list = zip(mT[::2], nexp[::2])
zipped_list = zip(mT, nexp)
sorted_list = sorted(zipped_list)
# sorted_list = sorted(zipped_list, reverse=True)
mT, nexp = np.array(sorted_list).T
# fin reordeno - ---------------
# parametrps -------------------


# reordeno - -------------------
# zipped_list = zip(mT[::2], nexp[::2])
zipped_list = zip(mT, nexp)
sorted_list = sorted(zipped_list)
# sorted_list = sorted(zipped_list, reverse=True)
mT, nexp = np.array(sorted_list).T


# cmap = matplotlib.cm.Wistia
# cmap = matplotlib.cm.autumn
# cmap = matplotlib.cm.autumn_r
cmap = matplotlib.cm.YlOrBr_r
# nexp = [24]

filename = f"0V_semiwidth_{semiancho}"
rango = (5, -8)

modulo = False
##############################################################################

# determino numero de filas y columnas del grafico-----------------------------

raiz = np.sqrt(len(nexp))
if int(raiz)**2 == len(nexp):
    filas = int(raiz)
    columnas = int(raiz)
else:
    filas = int(raiz//1)  # parte entera
    columnas = int(raiz//1)+1

# otro caso limitado
if filas*columnas < len(nexp):
    columnas += 1


size = 10
figura = plt.figure(1, figsize=(size, size*filas/columnas),
                    constrained_layout=False)
gs = figura.add_gridspec(filas, columnas)
axes = gs.subplots(sharex=True, sharey=True)
# ------------------------------------------------------------------------------

# chequeo que el semiancho no provoque solapamiento:---------------------------
ps = np.array(peaks)
dists = np.abs(ps - np.roll(ps, -1))
if any(dists < 2*semiancho):
    msj1 = f"El semiancho {semiancho} ppm hace que se solapen algunos peaks: {peaks} ppm"
    msj2 = f"\t El maximo semiancho posible es {np.min(dists/2)}"
    msj = f"{msj1}\n{msj2}"
    raise Exception(msj)
# ------------------------------------------------------------------------------
fig_specs, ax_specs = plt.subplots(nrows=len(peaks), ncols=len(peaks), num=78945975413)

espectros1d = []
integrales = []
espectros = []
for n in range(len(nexp)):
    directorio = f"{path}{int(nexp[n])}/"
    print(directorio, "  mT=", str(mT[n]))
    datos = DatosProcesados2D(directorio)
    datos.espectro.ppmSelect2D(rango)
    # print(datos.title)

    if modulo:
        re = datos.espectro.real
        im = datos.espectro.imag
        spec = np.abs(re+1j*im)
    else:
        spec = datos.espectro.real
        im = datos.espectro.imag

    ppm_x = datos.espectro.ppmAxis
    ppm_y = datos.espectro.ppmAxisInd
    ppm_x_grid = datos.espectro.ppmGridDir
    ppm_y_grid = datos.espectro.ppmGridInd

    vmax = np.max(spec)
    ax = axes[np.unravel_index(n, (filas, columnas))]
    ax.set_aspect('equal', 'box')
    print(f"mixing time: {mT[n]} ms")
    ax.set_title(f"mixing time: {mT[n]} ms")
    # plt.contour(ppm_x_grid, ppm_y_grid, spec, 30)
    # ax.contour(ppm_x_grid, ppm_y_grid, spec, 30, colors=['k'])
    ax.pcolormesh(ppm_x_grid, ppm_y_grid, spec, vmax=0.1*np.max(spec))

    ax.set_xlim((max(rango), min(rango)))
    ax.set_ylim((max(rango), min(rango)))

    Int_n = np.zeros([len(peaks), len(peaks)])
    for i in range(len(peaks)):
        for j in range(len(peaks)):
            xc = peaks[i]
            yc = peaks[j]

            ppmx_tmp = ppm_x[np.abs(ppm_x-xc) < semiancho]
            ppmy_tmp = ppm_y[np.abs(ppm_y-yc) < semiancho]

            xi = np.where(ppm_x == ppmx_tmp[0])[0][0]
            xf = np.where(ppm_x == ppmx_tmp[-1])[0][0]
            yi = np.where(ppm_y == ppmy_tmp[0])[0][0]
            yf = np.where(ppm_y == ppmy_tmp[-1])[0][0]

            spec_tmp = spec[yi:yf+1, xi:xf+1]

            #ax.contourf(ppmx_tmp, ppmy_tmp, spec_tmp, 5, cmap='inferno')            
            rect = Rectangle((xc-semiancho, yc-semiancho),
                              2*semiancho, 2*semiancho,
                              facecolor="none", edgecolor='r', 
                              linewidth=1)
            ax.add_patch(rect)
            ax.set_xticks([])
            ax.set_yticks([])

            Int_n[i, j] = Integrar(spec_tmp, x=ppmx_tmp, y=ppmy_tmp)

            # Grafico espectros 1D
            color_index = n / (len(mT) - 1) * 0.9  # Normalized index
            color = cmap(color_index)

            axx = ax_specs[i, j]
            if i == j:
                idx_y = find_nearest(ppm_y, peaks[i])
                spec1d = spec[idx_y, :]
                ppmaxis = ppm_y
                axx.set_title(f"Horizontal, pico {i+1}")

                # if i == 0:
                #     spec1dim = im[idx_y, :]
                #     # s1d, _ = autophase(spec1d+1j*spec1dim, x=ppm_x)
                #     #fase = fases[n]
                #     fase = 0
                #     s1d = (spec1d+1j*spec1dim)*np.exp(-1j*fase*np.pi/180)
                #     plt.figure(45138745953548)
                #     plt.plot(ppm_x, s1d,
                #              label=f"{mT[n]:.0f} ms", color=color)
                #     plt.xlim([max(rango), min(rango)])
                #     # plt.ylim([-200, 2000])
                #     plt.xlabel(r"$\Delta\delta$ [ppm]")
            else:
                idx_x = find_nearest(ppm_x, peaks[j])
                spec1d = spec[:, idx_x]
                ppmaxis = ppm_x
                axx.set_title(f"Vertical, pico {j+1}")
            axx.plot(spec1d, label=f"{mT[n]:.0f} ms", color=color)
            if j == 0 and i == 2:
                axx.legend(ncol=3)

    integrales.append(Int_n)
    # if mT[n]>=120: stop

    # %%
    # extraigo espectro el primer espectro de la exsy
    datos.set_fid()
    Fid1D = datos.fid
    Fid1D.signal = Fid1D.signal[0, :]  # selecciono la primera fid
    Fid1D.update_re_im()

    Espectro1D = Espectro(Fid1D)

    # espectro 1D como el primer espectro de la exsy
    espectros1d.append(Espectro1D.real)
    # espectros1d.append(np.sum(np.real(spec), axis=0)) ### espectro 1D como integral del espectro 2D

    # %%

integrales = np.array(integrales)
# plt.figure(101010)
# plt.subplot(filas, columnas,n+1)
# plt.xticks([])
# plt.yticks([])
# plt.contour(spec,30, cmap='inferno')
# plt.title("mixing time: {} ms".format(mT[n]))
# if guardar:
#     archivo_out = "mixTime_"+str(mT[n])+"ms"
#     np.savetxt(savepath+archivo_out+'.dat', datos.espectro.real)
#     np.savetxt(savepath+archivo_out+'_ppmDir.dat', datos.espectro.ppmAxis)
#     np.savetxt(savepath+archivo_out+'_ppmInd.dat', datos.espectro.ppmAxisInd)

# espectros.append(datos.espectro.real)
# X = exsy.ppmGridDir
# Y = exsy.ppmGridInd

# %%

for i in range(len(peaks)):
    for j in range(len(peaks)):
        if i < j:
            pii = integrales[:, i, i]
            pij = integrales[:, i, j]
            pji = integrales[:, j, i]
            pjj = integrales[:, j, j]
            # pii=1
            # pjj=1

            # cross peaks/ diag. peaks: medio, up y down
            cpu = pij/(pii+pjj)
            cpd = pji/(pii+pjj)
            cp = (pij+pji)/(pii+pjj)/2

            # CROSS PEAKS
            ii = i+1
            jj = j+1
            plt.figure(ii*10+jj)
            plt.title(f"peaks {ii} ({peaks[i]} ppm) and {jj} ({peaks[j]} ppm)")
            plt.plot(mT, cpu, 'ro-', label=f"P{ii}{jj}/(P{ii}+P{jj})")
            plt.plot(mT, cpd, 'bo-', label=f"P{jj}{ii}/(P{ii}+P{jj})")
            plt.plot(mT, cp, 'ko-',
                     label=f"[(P{ii}{jj}+P{jj}{ii})/2] / (P{ii}+P{jj})")
            # plt.plot(mT, (pii+pjj), 'ro-', label= f"(P{ii}+P{jj})")
            plt.xlabel("Mixing time [ms]")
            plt.ylabel("Cross Peaks/Diag. Peaks [a.u]")
            plt.legend()

            # plt.figure(ii*100+jj)
            # plt.title(f"peaks {ii} ({peaks[i]} ppm) and {jj} ({peaks[j]} ppm)")
            # plt.plot(mT, cp/cp[-1], 'o-', label= f"semiancho={semiancho}")
            # # plt.plot(mT, (pii+pjj), 'ro-', label= f"(P{ii}+P{jj})")
            # plt.xlabel("Mixing time [ms]")
            # plt.ylabel("Cross Peaks/Diag. Peaks [NORM]")
            # plt.legend()

            plt.figure(ii*1000+jj)
            plt.title(f"peaks {ii} ({peaks[i]} ppm) and {jj} ({peaks[j]} ppm)")
            plt.plot(mT, pij, 'ro-', label=f"UP")
            plt.plot(mT, pji, 'bo-', label=f"DOWN")
            # plt.plot(mT, (pii+pjj), 'ro-', label= f"(P{ii}+P{jj})")
            plt.xlabel("Mixing time [ms]")
            plt.ylabel("Cross Peaks")
            plt.legend()

            # guardado
            header = f"# peaks P{ii}:{peaks[i]} ppm and P{jj} {peaks[j]} ppm\n"
            # Up and Down son los que aparecen por arriba y por debajo de la diagonal
            header += "# CrossPeaks/Diag.Peaks     : Pomedio de ambos peaks cruzados\n"
            header += "# CrossPeaksUp/Diag.Peaks   : peaks cruzados por arriba de la diagonal (Pij con i<j)\n"
            header += "# CrossPeaksDown/Diag.Peaks : peaks cruzados por abajo de la diagonal  (Pij con i>j)\n"
            header += f"# MixingTime_[ms]\tCrossPeaks/Diag.Peaks\tCrossPeaksUp/Diag.Peaks\tCrossPeaksDown/Diag.Peaks\tP{i}{j}(up)\tP{j}{i}(down)\tP{ii}\tP{jj}"
            data = np.array([mT, cp, cpu, cpd, pij, pji, pii, pjj]).T
            np.savetxt(f"{savepath}/{filename}_peaks{ii}{jj}.dat",
                       data, header=header)

# %%
# DIAGONAL PEAKS
for i in range(len(peaks)):
    pii = integrales[:, i, i]
    plt.figure(100)
    # plt.plot(mT, pii/pii[0], 'o-', label= f"P{i+1}: ({peaks[i]} ppm)")
    plt.plot(mT, pii, 'o-', label=f"P{i+1}: ({peaks[i]} ppm)")
    plt.xlabel("Mixing time [ms]")
    plt.ylabel("Norm. Diag. Peaks")
    plt.legend()
# %%
# DIAGONAL PEAKS norm
for i in range(len(peaks)):
    pii = integrales[:, i, i]
    plt.figure(101)
    plt.plot(mT, pii/pii[0], 'o-', label=f"P{i+1}: ({peaks[i]} ppm)")
    plt.xlabel("Mixing time [ms]")
    plt.ylabel("Norm. Diag. Peaks")
    plt.legend()
# %% espectros 1d
# plt.figure(951623)
# ppm_x = Espectro1D.ppmAxis
# for n in range(len(espectros1d)):
#     spec1d = espectros1d[n]
#     # x_norm =
#     # norm_factor = spec1d[np.abs(ppm_x-x_norm)]
#     plt.plot(ppm_x, np.abs(spec1d)/np.max(np.abs(spec1d)),
#              label=f"mixing time: {mT[n]:.0f} ms")
# plt.legend()

# %% guardo matrices intensidad
np.savetxt(f"{savepath}/mixingTimes.dat", mT.T, header="mixig times [ms]")
for n in range(len(integrales)):
    matriz = integrales[n]
    filename = f'matrizI_mT{mT[n]:.0f}ms'
    data = np.array(matriz)
    np.savetxt(f"{savepath}/{filename}.dat", data)


# np.savetxt(f"{savepath}DATA2D_mT120.dat", spec)
# np.savetxt(f"{savepath}DATA2D_mT120_ppmAxis.dat", ppm_x)
# np.savetxt(f"{savepath}DATA2D_mT120_ppmAxisInd.dat", ppm_y)

# re = datos.espectro.real
# im = datos.espectro.imag
# mod = np.abs(re+1j*im)
# np.savetxt(f"{savepath}DATA2D_mT120_ABS.dat", mod)
# np.savetxt(f"{savepath}DATA2D_mT120_ABS.dat", mod)
