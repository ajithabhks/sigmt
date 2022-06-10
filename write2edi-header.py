# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 10:37:43 2021

@author: AJITHABH
"""
import pandas as pd

lat = loc.get('lat')
lat_d = np.floor(lat)
lat_m2 = (lat - lat_d) * 60
lat_m = np.floor(lat_m2)
lat_s = lat_m2 - lat_m
del lat_m2
lat_s = lat_s * 60
lon = loc.get('lon')
lon_d = np.floor(lon)
lon_m2 = (lon - lon_d) * 60
lon_m = np.floor(lon_m2)
lon_s = lon_m2 - lon_m
del lon_m2
lon_s = lon_s * 60

from datetime import date
today = date.today()
d1 = today.strftime("%m/%d/%Y")

f = open("C:/Users/Ajithabh/Desktop/Outputs/KL33A/KL33A-HEADER.edi", "x")
#===== Header =====
f.write(">HEAD\n")
f.write("  DATAID=" + procinfo.get('selectedsite'))
f.write("\n  ACQBY=" + "\"\"")
f.write("\n  FILEBY=" + "\"\"")
f.write("\n  ACQDATE=" + pd.to_datetime(timeline[0]-719529,unit='D').strftime("%m/%d/%Y"))
f.write("\n  FILEDATE=" + d1)
f.write("\n  PROSPECT=" + procinfo.get('selectedsite'))
f.write("\n  LAT=" + str(round(lat_d)) + ":" + str(round(lat_m)) + ":" + str(round(lat_s,2)))
f.write("\n  LONG=" + str(round(lon_d)) + ":" + str(round(lon_m)) + ":" + str(round(lon_s,2)))
f.write("\n  ELEV=" + str(loc.get('elev')))
f.write("\n  STDVERS= \"SEG 1.0\"")
f.write("\n  PROGVERS= ")
f.write("\n  PROGDATE= ")
f.write("\n  MAXSECT=999")
f.write("\n  EMPTY=1.0E32")

#===== INFO section =====
f.write("\n\n>INFO\n")
f.write("  MAXINFO=500")
f.write("\n\n")
f.write("------------------------------------------------------------")

#===== DEFINEMEAS section =====
f.write("\n\n>=DEFINEMEAS\n")
f.write("  MAXCHAN=" + "9")
f.write("\n  MAXRUN=" + "999")
f.write("\n  MAXMEAS=" + "1000")
f.write("\n  REFTYPE=" + "CART")
f.write("\n  REFLOC=" + procinfo.get('selectedsite'))
f.write("\n  REFLAT=" + str(round(lat_d)) + ":" + str(round(lat_m)) + ":" + str(round(lat_s,2)))
f.write("\n  REFLONG=" + str(round(lon_d)) + ":" + str(round(lon_m)) + ":" + str(round(lon_s,2)))
f.write("\n  REFELEV="+ str(loc.get('elev')))

#===== EMEAS & HMEAS section =====
f.write("\n\n>EMEAS  ID="+str(measid.get('Ex'))+"  CHTYPE=EX  X=-4.000000000E+01")
f.write("\n  Y=0.000000000E+00  Z=0.000000000E+00")
f.write("\n  ACQCHAN=\"ADU07/EFP06 /0/\"  GAIN=1  MEASDATE=" + pd.to_datetime(timeline[0]-719529,unit='D').strftime("%m/%d/%Y"))
f.write("\n  X2=4.000000000E+01  Y2=0.000000000E+00  Z2=0.000000000E+00")
f.write("\n>EMEAS  ID="+str(measid.get('Ey'))+"  CHTYPE=EY  X=0.000000000E+00")
f.write("\n  Y=0.000000000E+00  Z=0.000000000E+00")
f.write("\n  ACQCHAN=\"ADU07/EFP06 /0/\"  GAIN=1  MEASDATE=" + pd.to_datetime(timeline[0]-719529,unit='D').strftime("%m/%d/%Y"))
f.write("\n  X2=4.000000000E+01  Y2=0.000000000E+00  Z2=0.000000000E+00")

f.write("\n>HMEAS  ID="+str(measid.get('Hx'))+"  CHTYPE=HX  X=0.000000000E+00")
f.write("\n  Y=0.000000000E+00  Z=0.000000000E+00")
f.write("\n  ACQCHAN=\"ADU07/MFS06 /0/\"  GAIN=1  MEASDATE=" + pd.to_datetime(timeline[0]-719529,unit='D').strftime("%m/%d/%Y"))
f.write("\n  AZM=0.000000000E+00  DIP=0.000000000E+00")
f.write("\n  SENSOR=MFS06 /0")
f.write("\n>HMEAS  ID="+str(measid.get('Hy'))+"  CHTYPE=HY  X=0.000000000E+00")
f.write("\n  Y=0.000000000E+00  Z=0.000000000E+00")
f.write("\n  ACQCHAN=\"ADU07/MFS06 /0/\"  GAIN=1  MEASDATE=" + pd.to_datetime(timeline[0]-719529,unit='D').strftime("%m/%d/%Y"))
f.write("\n  AZM=9.000000000E+01  DIP=0.000000000E+00")
f.write("\n  SENSOR=MFS06 /0")
f.write("\n>HMEAS  ID="+str(measid.get('Hz'))+"  CHTYPE=HZ  X=0.000000000E+00")
f.write("\n  Y=0.000000000E+00  Z=0.000000000E+00")
f.write("\n  ACQCHAN=\"ADU07/MFS06 /0/\"  GAIN=1  MEASDATE=" + pd.to_datetime(timeline[0]-719529,unit='D').strftime("%m/%d/%Y"))
f.write("\n  AZM=0.000000000E+00  DIP=0.000000000E+00")
f.write("\n  SENSOR=MFS06 /0")

#===== MTSECT section =====
f.write("\n\n>=MTSECT\n")
f.write("  SECTID=" + procinfo.get('selectedsite'))
f.write("\n  NFREQ=" + str(np.size(cohEx)))
f.write("\n  HX=" + str(measid.get('Hx')))
f.write("\n  HY=" + str(measid.get('Hy')))
f.write("\n  HZ=" + str(measid.get('Hz')))
f.write("\n  EX=" + str(measid.get('Ex')))
f.write("\n  EY=" + str(measid.get('Ey')))
if 'measidR' in locals():
    f.write("\n  RX=" + str(measidR.get('Hx')))
    f.write("\n  RY=" + str(measidR.get('Hy')))


f.close()