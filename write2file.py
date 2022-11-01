# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 15:19:53 2020

@author: AJITHABH K. S.
Last modified: 21-07-2022

This script helps to create a text file with all impedance values, tipper data and variances.

This is required because, processing is done for data with different sampling frequencies. So,
text file need to be created with data of each target frequencies (For e.g. 1024 Hz). After
processing all data, join the text files to create a final text file contains data for all target
frequencies.

Then it can be converted to EDI format using 'write2edi-data.py' script.

Give path to save the text file with a filename in variable 'f'.

If the sampling frequency is 1024 Hz, give filename as '1024Hz.txt'

Read usermanual for more details.
"""

ZxxR = np.real(Z_huber.get('Zxx'))
ZxxI = np.imag(Z_huber.get('Zxx'))
ZxyR = np.real(Z_huber.get('Zxy'))
ZxyI = np.imag(Z_huber.get('Zxy'))
ZyxR = np.real(Z_huber.get('Zyx'))
ZyxI = np.imag(Z_huber.get('Zyx'))
ZyyR = np.real(Z_huber.get('Zyy'))
ZyyI = np.imag(Z_huber.get('Zyy'))
ZxxVar = Zvar.get('xx')
ZxyVar = Zvar.get('xy')
ZyxVar = Zvar.get('yx')
ZyyVar = Zvar.get('yy')
TxR = np.real(Tx)
TxI = np.imag(Tx)
TyR = np.real(Ty)
TyI = np.imag(Ty)
f = open('D:/NGRI/PROCESSING FOLDER/AMT/C6DAYPOFF/bands/65536Hz.txt', 'w')
for i in range(ftlist.shape[0]):
    f.write("%.9E \t %.9E  \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \t %.9E \n" % (ftlist[i][0],
    ZxxR[i][0], ZxxI[i][0], ZxyR[i][0], ZxyI[i][0],ZyxR[i][0], ZyxI[i][0],ZyyR[i][0], 
    ZyyI[i][0], ZxxVar[i][0], ZxyVar[i][0], ZyxVar[i][0], ZyyVar[i][0],
    TxR[i][0],TxI[i][0],TyR[i][0],TyI[i][0],TxVar[i][0],TyVar[i][0]))
f.close()
