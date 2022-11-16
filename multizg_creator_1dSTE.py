# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 18:33:41 2021

@author: santi
"""

import numpy as np
import shutil
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove

def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    # cuento cuantos reemplazos ocurrieron con count.
    count = 0
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_line = line.replace(pattern, subst)
                new_file.write(new_line)
                if new_line != line:
                  count += 1
                
    #Copy the file permissions from the old file to the new file
    copymode(file_path, abs_path)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)
    
    if count==0:
      warning = f"\nNo se encontro el patron:\n\t '{pattern}'\n"\
                f"En el archivo acqus de la capeta:\n\t {file_path}"
      raise Warning(warning)
    return count

# directorio de datos
path = "S:/Doctorado/Carbones/300MHz/2021-12-22_tmp/"
# directorio de guradado
# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/N32_Delta1.5ms/"

path = "S:/NMRdata/2022_Glicerol-agua_CNEA/2022-11-14_Diff_Silica_Agua-Glicerol-LiCl/"
       
expnum_i = 30
Nexp = 16 # numero total incluido la carpeta original


### crear nuevas carpetas
#expnum = np.arange(100,121)
for j in range(1,Nexp):
  try:
    shutil.copytree(path+str(expnum_i), path+str(expnum_i+j))
  except FileExistsError:
    pass
#%%
expnums = np.arange(expnum_i, expnum_i+Nexp)
# el gmin en este caso es 1, si quiero otro debo modificar esto
# glist = np.linspace(1, gmax, Nexp)
# glist = np.array([ 5.  , 22.7 , 31.71, 38.68, 44.57, 49.77, 54.47, 58.8 , 62.83,
#        66.62, 70.21, 73.62, 76.87, 80.  ])
# glist  = np.array([28.45,34.67,20.43,39.93,31.71,14.87,42.31,24.77,37.39])
# d1list = np.array([ 1.00, 1.00, 1.00, 1.30, 1.20, 1.00, 1.50, 1.00, 1.00])



glist = np.array([5.0, 43.18, 35.38, 25.26, 49.78, 55.6, 60.86, 65.71, 18.21, 39.47, 30.74, 58.29, 63.33, 46.6, 52.77, 68.0])
d1list = np.array([1.00, 2.00, 1.00, 1.50, 1.20, 1.50, 1.50, 2.00, 2.00, 1.00, 2.00, 2.00, 1.00, 1.30, 1.30, 1.50])

for n in range(glist.size):
  expn = expnums[n]
  g = glist[n]
  d1 = d1list[n]
  print(f"expn={expn}, g={g:.2f}, d1={d1}")
  file_path = path+'{}/acqu'.format(expn)
  
  # cambio el GPZ6
  busca     = "0 0 0 0 0 0 5 -17.13"
  reemplaza = "0 0 0 0 0 0 {:.2f} -17.13".format(g)
  replace(file_path, busca, reemplaza)
  
  # cambio el D1
  # busca     = "3e-06 1 0.02 3e-06 2e-05 0 1.5e-05"
  # reemplaza = "3e-06 {:.2f} 0.02 3e-06 2e-05 0 1.5e-05".format(d1)
  # replace(file_path, busca, reemplaza)

  # cambio el titulo
  file_path = path+'{}/pdata/1/title'.format(expn)
  busca     = 'gp=5.00%'
  reemplaza = 'gp={:.2f}%'.format(g)
  replace(file_path, busca, reemplaza)

#%%

# np.savetxt(path+"gp_list.dat", np.array([expnums, glist]).T, fmt="%i\t%.2f" ,header="expnum\tg (%)")