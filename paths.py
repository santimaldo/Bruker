# -*- coding: utf-8 -*-
"""
Rutas raiz de datos y de analisis, para no escribirlas a mano en cada script.

Uso en un script:
    from paths import NMRDATA, PROJECTS, NMRDATA_LOCAL
    path = f"{NMRDATA}/400dnp/2025-07-01_muestra/"

Valores por defecto (se calculan solos en cualquier PC con OneDrive):
    NMRDATA       -> <OneDrive>/NMRdata  (si no existe: NMRDATA_LOCAL)
    PROJECTS      -> <OneDrive>/Projects
    NMRDATA_LOCAL -> <home>/Documents/NMRdata

Para cambiarlos en una PC en particular, crear 'paths_local.json' en esta
carpeta (no se sube a git, ver paths_local.example.json), o definir las
variables de entorno BRUKER_NMRDATA, BRUKER_PROJECTS, BRUKER_NMRDATA_LOCAL.
Prioridad: variable de entorno > paths_local.json > valor por defecto.
"""
import json
import os
from pathlib import Path

_home = Path.home()
_onedrive = Path(os.environ.get("OneDriveCommercial")
                 or os.environ.get("OneDrive")
                 or _home / "OneDrive - University of Cambridge")

_rutas = {
    "NMRDATA": _onedrive / "NMRdata",
    "PROJECTS": _onedrive / "Projects",
    "NMRDATA_LOCAL": _home / "Documents" / "NMRdata",
}
# Los datos se mudaron de OneDrive a Documents: si no existe la carpeta de
# OneDrive, NMRDATA apunta a la copia local.
if not _rutas["NMRDATA"].exists():
    _rutas["NMRDATA"] = _rutas["NMRDATA_LOCAL"]

_local = Path(__file__).resolve().parent / "paths_local.json"
if _local.exists():
    with open(_local, "r") as f:
        _rutas.update(json.load(f))

for _nombre in _rutas:
    _rutas[_nombre] = str(os.environ.get(f"BRUKER_{_nombre}", _rutas[_nombre]))

NMRDATA = _rutas["NMRDATA"]
PROJECTS = _rutas["PROJECTS"]
NMRDATA_LOCAL = _rutas["NMRDATA_LOCAL"]


if __name__ == "__main__":
    for _nombre, _ruta in _rutas.items():
        print(f"{_nombre:14s} {_ruta}  {'(ok)' if Path(_ruta).exists() else '(NO EXISTE)'}")
