# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 18:33:41 2021

@author: santi
"""

import numpy as np
import shutil

# directorio de datos
path = "S:/Doctorado/LiMetal/116MHz/2021-04-05_testBB/"
# directorio de guradado
# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/N32_Delta1.5ms/"

expnum = np.arange(100,121)


for expn in expnum:
  shutil.copytree(path+str(expn), path+str(expn+100))


##### asi se reemplaza
with open('file.txt', 'r+') as f:
    
    #read file
    file_source = f.read()
    
    #replace 'PHP' with 'PYTHON' in the file
    replace_string = file_source.replace('PHP', 'PYTHON')
    
    #save output
    f.write(replace_string)