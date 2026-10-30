# -*- coding: utf-8 -*-
"""
Root folders used by every script. On another computer, edit these two lines.

Usage in a script::

    from paths import NMRDATA, HOME
    path = rf"{NMRDATA}\400dnp\2026-09-10_sample"
    savepath = rf"{HOME}\LiMetal\analysis"
"""

#: Bruker NMR data (spectrometer folders such as 400dnp, 500, 300neo).
NMRDATA = r"C:\Users\Santi\Documents\NMRdata"

#: Root of the analysis output folders.
HOME = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects"
