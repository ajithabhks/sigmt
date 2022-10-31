# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 15:47:38 2021

@author: AJITHABH K. S.
Last modified: 21-07-2022

This script can be used to arrange the impedance values, tipper and variances
in text file to as in EDI file.

See user manual for more detailed explanation

Give path of text file in 'txtfilename' variable.

Give path and EDI name to be created in the variable 'f'.

"""

import pandas as pd
import numpy as np

txtfilename = 'D:/NGRI/PROCESSING FOLDER/AMT/B4NIGHT/bands/all.txt'
f = open("D:/NGRI/PROCESSING FOLDER/AMT/B4NIGHT/B4NIGHT-DATA.edi", "x")

data = pd.read_csv(txtfilename, sep='\t', lineterminator='\n')
data = np.asarray(data)
freq = data[:,0]
zxxr = data[:,1]
zxxi = data[:,2]
zxyr = data[:,3]
zxyi = data[:,4]
zyxr = data[:,5]
zyxi = data[:,6]
zyyr = data[:,7]
zyyi = data[:,8]
zxxvar = data[:,9]
zxyvar = data[:,10]
zyxvar = data[:,11]
zyyvar  = data[:,12]
TxR = data[:,13]
TxI = data[:,14]
TyR = data[:,15]
TyI = data[:,16]
TxVar = data[:,17]
TyVar = data[:,18]

#===== FREQ section =====
f.write(">FREQ  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (freq[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== ZXXR section =====
f.write("\n \n")
f.write(">ZXXR  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zxxr[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== ZXXI section =====
f.write("\n \n")
f.write(">ZXXI  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zxxi[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== ZXX.VAR section =====
f.write("\n \n")
f.write(">ZXX.VAR  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zxxvar[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")

#===== ZXYR section =====
f.write("\n \n")
f.write(">ZXYR  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zxyr[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== ZXYI section =====
f.write("\n \n")
f.write(">ZXYI  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zxyi[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== ZXYI.VAR section =====
f.write("\n \n")
f.write(">ZXY.VAR  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zxyvar[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== ZYXR section =====
f.write("\n \n")
f.write(">ZYXR  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E" % (zyxr[i]))
    f.write("  ")
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== ZYXI section =====
f.write("\n \n")
f.write(">ZYXI  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zyxi[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")


#===== ZYX.VAR section =====
f.write("\n \n")
f.write(">ZYX.VAR  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zyxvar[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== ZYYR section =====
f.write("\n \n")
f.write(">ZYYR  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zyyr[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== ZYYI section =====
f.write("\n \n")
f.write(">ZYYI  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zyyi[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")


#===== ZYY.VAR section =====
f.write("\n \n")
f.write(">ZYY.VAR  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (zyyvar[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== TXR.EXP section =====
f.write("\n \n")
f.write(">TXR.EXP  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (TxR[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== TXI.EXP section =====
f.write("\n \n")
f.write(">TXI.EXP  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (TxI[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== TXVAR.EXP section =====
f.write("\n \n")
f.write(">TXVAR.EXP  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (TxVar[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== TYR.EXP section =====
f.write("\n \n")
f.write(">TYR.EXP  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (TyR[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== TYI.EXP section =====
f.write("\n \n")
f.write(">TYI.EXP  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (TyI[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")
        
#===== TYVAR.EXP section =====
f.write("\n \n")
f.write(">TYVAR.EXP  //"+str(np.size(freq)))
f.write("\n")

for i in range(np.size(freq)):
    k = i + 1
    f.write("%.9E " % (TyVar[i]))
    if (k%5 == 0) and (k != 0):
        f.write("\n")

f.write("\n \n")
f.write(">END")

f.close()