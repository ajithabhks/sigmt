# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 10:37:43 2021

@author: ajithabh

This script can be used to create EDI header information.

See user manual for more detailed explanation.

Give path and EDI name with 'HEADER' text in the variable 'f'. Then a header file
will be created with the information of the present site.

Join it with 'write2edi-data.py' output to get a final EDI file.

"""
from datetime import date
import pandas as pd

lat = loc.get('lat')
if lat > 0:
    lat_d = np.floor(lat)
    lat_m2 = (lat - lat_d) * 60
    lat_m = np.floor(lat_m2)
    lat_s = lat_m2 - lat_m
    del lat_m2
    lat_s = lat_s * 60
else:
    lat = -1 * lat
    lat_d = np.floor(lat)
    lat_m2 = (lat - lat_d) * 60
    lat_m = np.floor(lat_m2)
    lat_s = lat_m2 - lat_m
    del lat_m2
    lat_s = lat_s * 60
    lat_s = -1 * lat_s
lon = loc.get('lon')
if lon > 0:
    lon_d = np.floor(lon)
    lon_m2 = (lon - lon_d) * 60
    lon_m = np.floor(lon_m2)
    lon_s = lon_m2 - lon_m
    del lon_m2
    lon_s = lon_s * 60
else:
    lon = -1 * lon
    lon_d = np.floor(lon)
    lon_m2 = (lon - lon_d) * 60
    lon_m = np.floor(lon_m2)
    lon_s = lon_m2 - lon_m
    del lon_m2
    lon_s = lon_s * 60
    lon_d = -1 * lon_d

today = date.today()
d1 = today.strftime("%m/%d/%Y")

f = open("D:/NGRI/PROCESSING FOLDER/KL/KL33A/KL33A-HEADER.edi", "x")
# ===== Header =====
f.write(">HEAD\n")
f.write("  DATAID=" + procinfo.get('selectedsite'))
f.write("\n  ACQBY=" + "\"\"")
f.write("\n  FILEBY=" + "\"SigMT\"")
f.write("\n  ACQDATE=" +
        pd.to_datetime(timeline[0]-719529, unit='D').strftime("%m/%d/%Y"))
f.write("\n  FILEDATE=" + d1)
f.write("\n  PROSPECT=" + procinfo.get('selectedsite'))
f.write("\n  LAT=" + str(round(lat_d)) + ":" +
        str(round(lat_m)) + ":" + str(round(lat_s, 2)))
f.write("\n  LONG=" + str(round(lon_d)) + ":" +
        str(round(lon_m)) + ":" + str(round(lon_s, 2)))
f.write("\n  ELEV=" + str(loc.get('elev')))
f.write("\n  STDVERS= \"SEG 1.0\"")
f.write("\n  PROGVERS= v1.1.0")
f.write("\n  PROGDATE= ")
f.write("\n  MAXSECT=999")
f.write("\n  EMPTY=1.0E32")

# ===== INFO section =====
f.write("\n\n>INFO\n")
f.write("  MAXINFO=500")
f.write("\n\n")
f.write("------------------------------------------------------------")

# ===== DEFINEMEAS section =====
f.write("\n\n>=DEFINEMEAS\n")
f.write("  MAXCHAN=" + "9")
f.write("\n  MAXRUN=" + "999")
f.write("\n  MAXMEAS=" + "1000")
f.write("\n  REFTYPE=" + "CART")
f.write("\n  REFLOC=" + procinfo.get('selectedsite'))
f.write("\n  REFLAT=" + str(round(lat_d)) + ":" +
        str(round(lat_m)) + ":" + str(round(lat_s, 2)))
f.write("\n  REFLONG=" + str(round(lon_d)) + ":" +
        str(round(lon_m)) + ":" + str(round(lon_s, 2)))
f.write("\n  REFELEV=" + str(loc.get('elev')))

# ===== EMEAS & HMEAS section =====
f.write("\n\n>EMEAS ID="+"101.001" +
        " CHTYPE=EX X=0.000000E+00 Y=0.000000E+00 Z=0.000000E+00 X2=0.000000E+00 Y2=0.000000E+00 Z2=0.000000E+00")
f.write("\n>EMEAS ID="+"102.001" +
        " CHTYPE=EY X=0.000000E+00 Y=0.000000E+00 Z=0.000000E+00 X2=0.000000E+00 Y2=0.000000E+00 Z2=0.000000E+00")
f.write("\n>HMEAS ID="+"103.001" +
        " CHTYPE=HX X=0.000000E+00 Y=0.000000E+00 Z=0.000000E+00 X2=0.000000E+00 Y2=0.000000E+00 Z2=0.000000E+00")
f.write("\n>HMEAS ID="+"104.001" +
        " CHTYPE=HY X=0.000000E+00 Y=0.000000E+00 Z=0.000000E+00 X2=0.000000E+00 Y2=0.000000E+00 Z2=0.000000E+00")
f.write("\n>HMEAS ID="+"105.001" +
        " CHTYPE=HZ X=0.000000E+00 Y=0.000000E+00 Z=0.000000E+00 X2=0.000000E+00 Y2=0.000000E+00 Z2=0.000000E+00")
if 'measidR' in locals():
    f.write("\n>HMEAS ID="+"106.001" +
            " CHTYPE=HX X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 Y2=0.000000e+00 Z2=0.000000e+00")
    f.write("\n>HMEAS ID="+"107.001" +
            " CHTYPE=HY X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 Y2=0.000000e+00 Z2=0.000000e+00")

# ===== MTSECT section =====
f.write("\n\n>=MTSECT\n")
f.write("  SECTID=" + procinfo.get('selectedsite'))
f.write("\n  NFREQ=" + str(np.size(cohEx)))
f.write("\n  EX=" + "101.001")
f.write("\n  EY=" + "102.001")
f.write("\n  HX=" + "103.001")
f.write("\n  HY=" + "104.001")
f.write("\n  HZ=" + "105.001")
if 'measidR' in locals():
    f.write("\n  RX=" + "106.001")
    f.write("\n  RY=" + "107.001")


f.close()
