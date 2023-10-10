import nmrglue as ng
import matplotlib.pyplot as plt
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
        i = (integrate.simps(integrate.simps(Matriz)))
    else:
        i = (integrate.simps(integrate.simps(Matriz, x=x), x=y))
    return i


# ###############################################################################
# path  ="S:/Doctorado/Carbones/300MHz/2019-10-24_Carbones_MAS_EXSY-reanalisis2022/"
# savepath = "S:/temp/"
# # mT =   [1,1,5,10,25,50,75,100] #ms
# # nexp = [8,9,5, 7,15,10,13, 14]
# mT =   [1,25,50,75,100] #ms
# nexp = [9,15,10,13, 14]
# # fin reordeno - ---------------
# picos = [-0.5,-4.3]
# # parametrps -------------------
# # semiancho de integracion (ppm)
# semiancho = 1
# rango = (-12,6)
# ###############################################################################


##############################################################################
# CMK3ACTB Li
# path  ="S:/Doctorado/Carbones/300MHz/2022-05-19_Carbones_CMK3ACT/"
# path = "S:/NMRdata/2021_Carbones_Sofi/2022-05-19_Carbones_CMK3ACT/"  # compu Ofi
# savepath = path
# filename = "7Li_EXSY"
# mT = [1000, 300, 100, 600, 1, 25, 1250, 50,
#       450, 200, 800, 10, 150, 1500, 700, 375]
# nexp = np.arange(5, 5+2*len(mT), 2)
# # reordeno - -------------------
# zipped_list = zip(mT, nexp)
# sorted_list = sorted(zipped_list)
# mT, nexp = np.array(sorted_list).T
# # fin reordeno - ---------------
# picos = [-1.57, -3.75]
# # picos   = [-1.6,-2.8,-3.7]
# # parametrps -------------------
# # semiancho de integracion (ppm)
# semiancho = 0.5
# rango = (-7, 0)
# modulo = False
###############################################################################

# #############################################################################
# #############################################################################
# ###################CMK3ACTB 1H  # muestra "B". con cantidad de bulk correcta
path = "S:/Doctorado/Carbones/300MHz/2022-05-12_Carbones_CMK3ACT/"  # Acer
# path = "S:/NMRdata/2021_Carbones_Sofi/2022-05-12_Carbones_CMK3ACT/"  # compu Ofi
# savepath = "S:/tmp"
savepath = "S:/tmp/"
mT = [1, 100, 350, 1000, 10, 35, 600, 5, 20, 200, 75, 275,
      50, 800, 150, 700, 900, 500, 420, 120, 1200, 1500]  # 1H
nexp = np.arange(22, 22+2*len(mT), 2)
# reordeno - -------------------
zipped_list = zip(mT[::2], nexp[::2])
# zipped_list = zip(mT, nexp)
sorted_list = sorted(zipped_list)
# sorted_list = sorted(zipped_list, reverse=True)
mT, nexp = np.array(sorted_list).T
# fin reordeno - ---------------
picos = [3.61, 0.63, -3.3]
# parametrps -------------------
# semiancho de integracion (ppm)
semiancho = 0.5

# cmap = matplotlib.cm.Wistia
# cmap = matplotlib.cm.autumn
# cmap = matplotlib.cm.autumn_r
cmap = matplotlib.cm.YlOrBr_r
# nexp = [24]

filename = f"1H_EXSY_CMK3-ACT_semiancho{semiancho}"
rango = (-10, 10)

modulo = False
##############################################################################


# #############################################################################
# ################### 1H CMK3 sin activar
# #ATENCIOOOON: NO IMPLEMENTADO AUN!!!!
# path = "S:/Doctorado/Carbones/300MHz/2022-05-12_Carbones_CMK3ACT/"  # Acer
# path = "S:/NMRdata/2021_Carbones_Sofi/2021-12-27_Carbones_MAS/"  # compu Ofi
# # savepath = "S:/tmp"
# savepath = "S:/tmp"
# mT = [200, 1000, 100, 50, 10, 500, 20]  # 1H
# nexp = np.arange(21, 28)
# # reordeno - -------------------
# zipped_list = zip(mT, nexp)
# sorted_list = sorted(zipped_list)
# mT, nexp = np.array(sorted_list).T
# # fin reordeno - ---------------
# picos = [3.61, 0.63, -3.3]
# # parametrps -------------------
# # semiancho de integracion (ppm)
# semiancho = 0.5
# filename = f"1H_EXSY_CMK3_semiancho{semiancho}"
# rango = (-10, 10)

# modulo = False
# #############################################################################

# ###CMK3ACTB 1H  ### muestra "A". con MUCHA cantidad de bulk
# path  ="S:/CarbonesSofi/300MHz/2022-05-06_Carbones_CMK3ACT/" # compu Ofi
# savepath = "S:/CarbonesSofi/Analisis/2022-09_EXSY/"

# mT = [100,10,1000,200,20,500,50,150,5,350,750,1,900,600,275,75,675,420,120,35]
# nexp = np.arange(31,31+len(mT))
# # reordeno - -------------------
# zipped_list = zip(mT, nexp)
# sorted_list = sorted(zipped_list)
# mT, nexp = np.array(sorted_list).T
# # fin reordeno - ---------------
# picos = [3.5, 0.6]
# # parametrps -------------------
# # semiancho de integracion (ppm)
# semiancho = 0.5


# filename = f"1H_EXSY_CMK3-ACT-A_MuchoBulk_semiancho{semiancho}"
# rango = (-5,5)

# modulo=False
# ##########################################

##########################################
# M4 carbones HOracio  carbones SIN AGUA
# path  = "S:/Doctorado/Carbones/300MHz/2019-10-24_Carbones_MAS_EXSY-reanalisis2022/"
# savepath = "S:/temp/"

# mT =   [1,  5, 10, 25, 50, 75, 100]
# nexp = [9,  5,  7, 15, 10, 13,  14]
# # reordeno - -------------------
# zipped_list = zip(mT, nexp)
# sorted_list = sorted(zipped_list)
# mT, nexp = np.array(sorted_list).T
# # fin reordeno - ---------------
# picos = [0,-4.3]
# # parametrps -------------------
# # semiancho de integracion (ppm)
# semiancho = 1

# filename = f"1H_EXSY_M4_semiancho{semiancho}"
# rango = (-8,2)

# modulo=False
# ###############################################################################


# ###############################################################################
# M4 carbones HOracio    CARBON Y agua
# path = "S:/Doctorado/Carbones/300MHz/2019-10-24_Carbones_MAS_EXSY-reanalisis2022/"
# savepath = "S:/temp/"

# mT = [1, 5, 10, 25, 50, 100]
# nexp = [25, 29, 26, 31, 27, 30]
# # reordeno - -------------------
# zipped_list = zip(mT, nexp)
# sorted_list = sorted(zipped_list)
# mT, nexp = np.array(sorted_list).T
# # fin reordeno - ---------------
# picos = [4.5, 1.8, -2.8]
# # parametrps -------------------
# # semiancho de integracion (ppm)
# semiancho = 1

# filename = f"1H_EXSY_M4_semiancho{semiancho}"
# rango = (-5, 8)

# modulo = False
# ###############################################################################

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
ps = np.array(picos)
dists = np.abs(ps - np.roll(ps, -1))
if any(dists < 2*semiancho):
    msj1 = f"El semiancho {semiancho} ppm hace que se solapen algunos picos: {picos} ppm"
    msj2 = f"\t El maximo semiancho posible es {np.min(dists/2)}"
    msj = f"{msj1}\n{msj2}"
    raise Exception(msj)
# ------------------------------------------------------------------------------
fig_specs, ax_specs = plt.subplots(nrows=3, ncols=3, num=78945975413)

espectros1d = []
integrales = []
espectros = []
for n in range(len(nexp)):

    print(path+str(nexp[n])+"/", "  mT=", str(mT[n]))
    directorio = f"{path}{nexp[n]}/"
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

    ax.set_xlim((rango[1], rango[0]))
    ax.set_ylim((rango[1], rango[0]))

    Int_n = np.zeros([len(picos), len(picos)])
    for i in range(len(picos)):
        for j in range(len(picos)):
            xc = picos[i]
            yc = picos[j]

            ppmx_tmp = ppm_x[np.abs(ppm_x-xc) < semiancho]
            ppmy_tmp = ppm_y[np.abs(ppm_y-yc) < semiancho]

            xi = np.where(ppm_x == ppmx_tmp[0])[0][0]
            xf = np.where(ppm_x == ppmx_tmp[-1])[0][0]
            yi = np.where(ppm_y == ppmy_tmp[0])[0][0]
            yf = np.where(ppm_y == ppmy_tmp[-1])[0][0]

            spec_tmp = spec[yi:yf+1, xi:xf+1]

            ax.contourf(ppmx_tmp, ppmy_tmp, spec_tmp, 5, cmap='inferno')
            ax.set_xticks([])
            ax.set_yticks([])

            Int_n[i, j] = Integrar(spec_tmp, x=ppmx_tmp, y=ppmy_tmp)

            # Grafico espectros 1D
            color_index = n / (len(mT) - 1) * 0.9  # Normalized index
            color = cmap(color_index)

            axx = ax_specs[i, j]
            if i == j:
                idx_y = find_nearest(ppm_y, picos[i])
                spec1d = spec[idx_y, :]
                ppmaxis = ppm_y
                axx.set_title(f"Horizontal, pico {i+1}")

                if i == 0:
                    spec1dim = im[idx_y, :]
                    s1d, _ = autophase(spec1d+1j*spec1dim, x=ppm_x)
                    plt.figure(45138745953548)
                    plt.plot(ppm_x-3.6, s1d,
                             label=f"{mT[n]:.0f} ms", color=color)
                    plt.xlim([5, -15])
                    plt.ylim([-200, 2000])
                    plt.xlabel(r"$\Delta\delta$ [ppm]")
            else:
                idx_x = find_nearest(ppm_x, picos[j])
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

for i in range(len(picos)):
    for j in range(len(picos)):
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
            plt.title(f"Picos {ii} ({picos[i]} ppm) y {jj} ({picos[j]} ppm)")
            plt.plot(mT, cpu, 'ro-', label=f"P{ii}{jj}/(P{ii}+P{jj})")
            plt.plot(mT, cpd, 'bo-', label=f"P{jj}{ii}/(P{ii}+P{jj})")
            plt.plot(mT, cp, 'ko-',
                     label=f"[(P{ii}{jj}+P{jj}{ii})/2] / (P{ii}+P{jj})")
            # plt.plot(mT, (pii+pjj), 'ro-', label= f"(P{ii}+P{jj})")
            plt.xlabel("Mixing time [ms]")
            plt.ylabel("Cross Peaks/Diag. Peaks [a.u]")
            plt.legend()

            # plt.figure(ii*100+jj)
            # plt.title(f"Picos {ii} ({picos[i]} ppm) y {jj} ({picos[j]} ppm)")
            # plt.plot(mT, cp/cp[-1], 'o-', label= f"semiancho={semiancho}")
            # # plt.plot(mT, (pii+pjj), 'ro-', label= f"(P{ii}+P{jj})")
            # plt.xlabel("Mixing time [ms]")
            # plt.ylabel("Cross Peaks/Diag. Peaks [NORM]")
            # plt.legend()

            plt.figure(ii*1000+jj)
            plt.title(f"Picos {ii} ({picos[i]} ppm) y {jj} ({picos[j]} ppm)")
            plt.plot(mT, pij, 'ro-', label=f"UP")
            plt.plot(mT, pji, 'bo-', label=f"DOWN")
            # plt.plot(mT, (pii+pjj), 'ro-', label= f"(P{ii}+P{jj})")
            plt.xlabel("Mixing time [ms]")
            plt.ylabel("Cross Peaks")
            plt.legend()

            # guardado
            header = f"# Picos P{ii}:{picos[i]} ppm y P{jj} {picos[j]} ppm\n"
            # Up y Down son los que aparecen por arriba y por debajo de la diagonal
            header += "# CrossPeaks/Diag.Peaks     : Pomedio de ambos picos cruzados\n"
            header += "# CrossPeaksUp/Diag.Peaks   : Picos cruzados por arriba de la diagonal (Pij con i<j)\n"
            header += "# CrossPeaksDown/Diag.Peaks : Picos cruzados por abajo de la diagonal  (Pij con i>j)\n"
            header += f"# MixingTime_[ms]\tCrossPeaks/Diag.Peaks\tCrossPeaksUp/Diag.Peaks\tCrossPeaksDown/Diag.Peaks\tP{i}{j}(up)\tP{j}{i}(down)\tP{ii}\tP{jj}"
            data = np.array([mT, cp, cpu, cpd, pij, pji, pii, pjj]).T
            np.savetxt(f"{savepath}/{filename}_Picos{ii}{jj}.dat",
                       data, header=header)

# %%
# DIAGONAL PEAKS
for i in range(len(picos)):
    pii = integrales[:, i, i]
    plt.figure(100)
    # plt.plot(mT, pii/pii[0], 'o-', label= f"P{i+1}: ({picos[i]} ppm)")
    plt.plot(mT, pii, 'o-', label=f"P{i+1}: ({picos[i]} ppm)")
    plt.xlabel("Mixing time [ms]")
    plt.ylabel("Norm. Diag. Peaks")
    plt.legend()
# %%
# DIAGONAL PEAKS norm
for i in range(len(picos)):
    pii = integrales[:, i, i]
    plt.figure(101)
    plt.plot(mT, pii/pii[0], 'o-', label=f"P{i+1}: ({picos[i]} ppm)")
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
