# -*- coding: utf-8 -*-
"""
Funciones compartidas por los scripts de Kpopt (barrido de campo con bsms)
y de popt (subcarpetas 99n de TopSpin).

Antes estaban copiadas en cada script; importarlas desde aca:
    from kpopt_utils import read_bsms_field, shift_spectrum_by_bsms
    from kpopt_utils import find_popt_subfolders, extract_popt_parameters
"""
import ast
import os
import re
import numpy as np


################## Kpopt / bsms ###########################
def read_bsms_field(path_archivo):
    """
    Lee el archivo Klog_opt de la carpeta 'path_archivo' y devuelve los
    valores de bsms (la linea que empieza con 'x=').

    Parametros
    ----------
    path_archivo : str
        Carpeta del experimento (con o sin '/' final).

    Retorna
    -------
    np.ndarray con los valores de bsms.
    """
    with open(os.path.join(path_archivo, 'Klog_opt'), 'r') as f:
        for linea in f:
            if linea.startswith("x="):
                contenido = linea.strip().split("=", 1)[1]
                return np.array(ast.literal_eval(contenido))
    raise ValueError("No se encontro una linea que empiece con 'x='.")


def shift_spectrum_by_bsms(spectrum, ppmAxis, bsms, shift_factor=-1.1665e-2):
    """
    Desplaza un espectro en puntos enteros de la discretización ppm.

    Parámetros
    ----------
    spectrum : np.array
        Espectro 1D (real o complejo).
    ppmAxis : np.array
        Eje ppm original.
    bsms : float
        Valor de bsms correspondiente al espectro.
    shift_factor : float
        Factor de conversión bsms -> ppm.

    Retorna
    -------
    shifted_spectrum : np.array
        Espectro desplazado en puntos enteros.
    """
    # resolución digital del eje ppm
    dppm = np.abs(ppmAxis[1] - ppmAxis[0])
    # corrimiento continuo esperado
    shift_ppm = shift_factor * bsms
    # cantidad entera de puntos digitales
    n_shift = int(np.round(shift_ppm / dppm))
    # desplazamiento sin condiciones periódicas
    shifted_spectrum = np.zeros_like(spectrum)
    if n_shift > 0:
        shifted_spectrum[n_shift:] = spectrum[:-n_shift]
    elif n_shift < 0:
        shifted_spectrum[:n_shift] = spectrum[-n_shift:]
    else:
        shifted_spectrum = spectrum.copy()
    return shifted_spectrum


################## popt (subcarpetas 99n) ###########################
def find_popt_subfolders(expn, path):
    """
    Find all subfolders that start with the given 
    experiment number followed by 99n:
    999, 998, ..., 1000-n.
    """
    experiment_numbers = []
    for folder in os.listdir(path):
        if folder.startswith(f"{expn}"):
            try:
                suffix = int(folder[len(f"{expn}"):])
                experiment_numbers.append(suffix)
            except ValueError:
                continue
    return np.array(sorted(experiment_numbers, reverse=True))


def extract_popt_parameters(subfile, parameter, path=".", expn=1):
    """
    Extracts the specified parameter from the popt.protocol file.
    The parameter is specified as a string.

    typical options are:
    "P 15" for CP contact time
    """
    filename = f"{path}/{expn}/popt.protocol.{subfile}"
    data = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        start_reading = False
        column_index = None
        for line in lines:
            if line.startswith("Experiment") and parameter in line:
                headers = re.split(r'\s{2,}', line)
                column_index = headers.index(parameter)
                start_reading = True
                continue
            if start_reading:
                columns = line.split()
                if columns == [] or columns[0].isdigit() == False:
                    break
                
                data.append(float(columns[column_index]))
    return np.array(data)
