import math
from datetime import date

import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QFileDialog

from sigmt.__version__ import __version__


def save_edi(estimates, procinfo, project_setup):
    """
    Docs
    """

    lat = procinfo['lat']
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
    lon = procinfo['lon']
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

    # Data
    freq = estimates.frequency.values
    # For MT Only and MT + Tipper cases
    if project_setup['processing_mode'] == 'MT + Tipper' or project_setup['processing_mode'] == 'MT Only':
        zxx_r = np.real(estimates['zxx'].values)
        zxx_i = np.imag(estimates['zxx'].values)
        zxy_r = np.real(estimates['zxy'].values)
        zxy_i = np.imag(estimates['zxy'].values)
        zyx_r = np.real(estimates['zyx'].values)
        zyx_i = np.imag(estimates['zyx'].values)
        zyy_r = np.real(estimates['zyy'].values)
        zyy_i = np.imag(estimates['zyy'].values)
        #
        zxx_var = estimates['zxx_var'].values
        zxy_var = estimates['zxy_var'].values
        zyx_var = estimates['zyx_var'].values
        zyy_var = estimates['zyy_var'].values
        #
        coh_ex = estimates['coh_ex'].values
        coh_ey = estimates['coh_ey'].values

    # For Tipper only and MT + Tipper cases
    if project_setup['processing_mode'] == 'MT + Tipper' or project_setup['processing_mode'] == 'Tipper only':
        tx_r = np.real(estimates['tzx'].values)
        tx_i = np.imag(estimates['tzx'].values)
        ty_r = np.real(estimates['tzy'].values)
        ty_i = np.imag(estimates['tzy'].values)
        tx_var = estimates['tzx_var'].values
        ty_var = estimates['tzy_var'].values

    today = date.today()
    today = today.strftime("%m/%d/%Y")
    #
    file_formats = "EDI Files (*.edi);;All Files (*)"
    edi_file, _ = QFileDialog.getSaveFileName(None, "Save File", "", file_formats)
    if edi_file:  # If user selects a path
        with open(edi_file, 'w', encoding='utf-8') as file_handle:
            file_handle.write(">HEAD\n")
            file_handle.write("  DATAID=" + "\"" + procinfo['localsite'] + "\"")
            file_handle.write("\n  ACQBY=" + "\"" + project_setup['acquired_by'] + "\"")
            file_handle.write("\n  FILEBY=" + "\"SigMT\"")
            file_handle.write("\n  ACQDATE=" + pd.to_datetime(procinfo['start_time'], unit='s').strftime("%m/%d/%Y"))
            file_handle.write("\n  FILEDATE=" + today)
            file_handle.write(f"\n  PROSPECT=\"{project_setup['project_name']}\"")
            file_handle.write(
                "\n  LAT=" + str(int(lat_d)) + ":" + str(int(np.floor(lat_m))) + ":" + str(round(lat_s, 2)))
            file_handle.write(
                "\n  LONG=" + str(int(lon_d)) + ":" + str(int(np.floor(lon_m))) + ":" + str(round(lon_s, 2)))
            file_handle.write("\n  ELEV=" + str(procinfo['elev']))
            file_handle.write("\n  STDVERS=\"SEG 1.0\"")
            file_handle.write("\n  PROGVERS=\"" + __version__ + "\"")
            file_handle.write("\n  PROGDATE= ")
            file_handle.write("\n  MAXSECT=999")
            file_handle.write("\n  EMPTY=1.000000000E+32")

            # ===== INFO section =====
            file_handle.write("\n\n>INFO\n")
            file_handle.write("  MAXINFO=500")
            file_handle.write("\n\n")
            file_handle.write("------------------------------------------------------------")

            # ===== DEFINEMEAS section =====
            file_handle.write("\n\n>=DEFINEMEAS\n")
            file_handle.write("  REFLOC=" + "\"" + procinfo['localsite'] + "\"")
            file_handle.write("\n  REFLAT=" + str(round(lat_d)) + ":" + str(round(lat_m)) + ":" + str(round(lat_s, 2)))
            file_handle.write("\n  REFLONG=" + str(round(lon_d)) + ":" + str(round(lon_m)) + ":" + str(round(lon_s, 2)))
            file_handle.write("\n  REFELEV=" + str(procinfo['elev']))

            if project_setup['processing_mode'] == 'MT + Tipper':
                file_handle.write(
                    "\n"
                    "\n>EMEAS ID=" + "101.001" + " CHTYPE=EX X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                file_handle.write(
                    "\n>EMEAS ID=" + "102.001" + " CHTYPE=EY X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                file_handle.write(
                    "\n>HMEAS ID=" + "103.001" + " CHTYPE=HX X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                file_handle.write(
                    "\n>HMEAS ID=" + "104.001" + " CHTYPE=HY X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                file_handle.write(
                    "\n>HMEAS ID=" + "105.001" + " CHTYPE=HZ X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                if procinfo['remotesite'] is not None:
                    file_handle.write(
                        "\n>HMEAS ID=" + "106.001" + " CHTYPE=HX X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                     "Y2=0.000000e+00 Z2=0.000000e+00")
                    file_handle.write(
                        "\n>HMEAS ID=" + "107.001" + " CHTYPE=HY X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                     "Y2=0.000000e+00 Z2=0.000000e+00")
            elif project_setup['processing_mode'] == 'MT only':
                file_handle.write(
                    "\n"
                    "\n>EMEAS ID=" + "101.001" + " CHTYPE=EX X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                file_handle.write(
                    "\n>EMEAS ID=" + "102.001" + " CHTYPE=EY X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                file_handle.write(
                    "\n>HMEAS ID=" + "103.001" + " CHTYPE=HX X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                file_handle.write(
                    "\n>HMEAS ID=" + "104.001" + " CHTYPE=HY X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                if procinfo['remotesite'] is not None:
                    file_handle.write(
                        "\n>HMEAS ID=" + "106.001" + " CHTYPE=HX X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                     "Y2=0.000000e+00 Z2=0.000000e+00")
                    file_handle.write(
                        "\n>HMEAS ID=" + "107.001" + " CHTYPE=HY X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                     "Y2=0.000000e+00 Z2=0.000000e+00")
            elif project_setup['processing_mode'] == 'Tipper only':
                file_handle.write(
                    "\n>HMEAS ID=" + "103.001" + " CHTYPE=HX X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                file_handle.write(
                    "\n>HMEAS ID=" + "104.001" + " CHTYPE=HY X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                 "Y2=0.000000e+00 Z2=0.000000e+00")
                if procinfo['remotesite'] is not None:
                    file_handle.write(
                        "\n>HMEAS ID=" + "106.001" + " CHTYPE=HX X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                     "Y2=0.000000e+00 Z2=0.000000e+00")
                    file_handle.write(
                        "\n>HMEAS ID=" + "107.001" + " CHTYPE=HY X=0.000000e+00 Y=0.000000e+00 Z=0.000000e+00 X2=0.000000e+00 "
                                                     "Y2=0.000000e+00 Z2=0.000000e+00")

            # ===== MTSECT section =====
            file_handle.write("\n\n>=MTSECT\n")
            file_handle.write("  SECTID=" + "\"" + procinfo['localsite'] + "\"")
            file_handle.write("\n  NFREQ=" + str(len(estimates.frequency)))
            if project_setup['processing_mode'] == 'MT + Tipper':
                file_handle.write("\n  EX=" + "101.001")
                file_handle.write("\n  EY=" + "102.001")
                file_handle.write("\n  HX=" + "103.001")
                file_handle.write("\n  HY=" + "104.001")
                file_handle.write("\n  HZ=" + "105.001")
                if procinfo['remotesite'] is not None:
                    file_handle.write("\n  RX=" + "106.001")
                    file_handle.write("\n  RY=" + "107.001")
            elif project_setup['processing_mode'] == 'MT Only':
                file_handle.write("\n  EX=" + "101.001")
                file_handle.write("\n  EY=" + "102.001")
                file_handle.write("\n  HX=" + "103.001")
                file_handle.write("\n  HY=" + "104.001")
                if procinfo['remotesite'] is not None:
                    file_handle.write("\n  RX=" + "106.001")
                    file_handle.write("\n  RY=" + "107.001")
            elif project_setup['processing_mode'] == 'Tipper only':
                file_handle.write("\n  HX=" + "103.001")
                file_handle.write("\n  HY=" + "104.001")
                file_handle.write("\n  HZ=" + "105.001")
                if procinfo['remotesite'] is not None:
                    file_handle.write("\n  RX=" + "106.001")
                    file_handle.write("\n  RY=" + "107.001")

            # ===== FREQ section =====
            file_handle.write("\n \n")
            file_handle.write(">FREQ  //" + str(np.size(freq)))
            file_handle.write("\n")

            for i in range(np.size(freq)):
                k = i + 1
                if freq[i] < 0:
                    file_handle.write("%.9E " % (freq[i]))
                else:
                    file_handle.write(" %.9E " % (freq[i]))
                if (k % 6 == 0) and (k != 0):
                    file_handle.write("\n")

            if project_setup['processing_mode'] == 'MT + Tipper' or project_setup['processing_mode'] == 'MT Only':
                # ===== ZXXR section =====
                file_handle.write("\n \n")
                file_handle.write(">ZXXR  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zxx_r[i] < 0:
                        file_handle.write("%.9E " % float(zxx_r[i]))
                    else:
                        file_handle.write(" %.9E " % float(zxx_r[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZXXI section =====
                file_handle.write("\n \n")
                file_handle.write(">ZXXI  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zxx_i[i] < 0:
                        file_handle.write("%.9E " % float(zxx_i[i]))
                    else:
                        file_handle.write(" %.9E " % float(zxx_i[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZXX.VAR section =====
                file_handle.write("\n \n")
                file_handle.write(">ZXX.VAR  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zxx_var[i] < 0:
                        file_handle.write("%.9E " % (zxx_var[i]))
                    else:
                        file_handle.write(" %.9E " % (zxx_var[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZXYR section =====
                file_handle.write("\n \n")
                file_handle.write(">ZXYR  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zxy_r[i] < 0:
                        file_handle.write("%.9E " % float(zxy_r[i]))
                    else:
                        file_handle.write(" %.9E " % float(zxy_r[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZXYI section =====
                file_handle.write("\n \n")
                file_handle.write(">ZXYI  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zxy_i[i] < 0:
                        file_handle.write("%.9E " % float(zxy_i[i]))
                    else:
                        file_handle.write(" %.9E " % float(zxy_i[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZXYI.VAR section =====
                file_handle.write("\n \n")
                file_handle.write(">ZXY.VAR  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zxy_var[i] < 0:
                        file_handle.write("%.9E " % (zxy_var[i]))
                    else:
                        file_handle.write(" %.9E " % (zxy_var[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZYXR section =====
                file_handle.write("\n \n")
                file_handle.write(">ZYXR  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zyx_r[i] < 0:
                        file_handle.write("%.9E " % float(zyx_r[i]))
                    else:
                        file_handle.write(" %.9E " % float(zyx_r[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZYXI section =====
                file_handle.write("\n \n")
                file_handle.write(">ZYXI  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zyx_i[i] < 0:
                        file_handle.write("%.9E " % float(zyx_i[i]))
                    else:
                        file_handle.write(" %.9E " % float(zyx_i[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZYX.VAR section =====
                file_handle.write("\n \n")
                file_handle.write(">ZYX.VAR  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zyx_var[i] < 0:
                        file_handle.write("%.9E " % float(zyx_var[i]))
                    else:
                        file_handle.write(" %.9E " % float(zyx_var[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZYYR section =====
                file_handle.write("\n \n")
                file_handle.write(">ZYYR  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zyy_r[i] < 0:
                        file_handle.write("%.9E " % float(zyy_r[i]))
                    else:
                        file_handle.write(" %.9E " % float(zyy_r[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZYYI section =====
                file_handle.write("\n \n")
                file_handle.write(">ZYYI  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zyy_i[i] < 0:
                        file_handle.write("%.9E " % float(zyy_i[i]))
                    else:
                        file_handle.write(" %.9E " % float(zyy_i[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== ZYY.VAR section =====
                file_handle.write("\n \n")
                file_handle.write(">ZYY.VAR  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if zyy_var[i] < 0:
                        file_handle.write("%.9E " % float(zyy_var[i]))
                    else:
                        file_handle.write(" %.9E " % float(zyy_var[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

            if project_setup['processing_mode'] == 'MT + Tipper' or project_setup['processing_mode'] == 'Tipper Only':
                # ===== TXR.EXP section =====
                file_handle.write("\n \n")
                file_handle.write(">TXR.EXP  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if tx_r[i] < 0:
                        file_handle.write("%.9E " % float(tx_r[i]))
                    else:
                        file_handle.write(" %.9E " % float(tx_r[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== TXI.EXP section =====
                file_handle.write("\n \n")
                file_handle.write(">TXI.EXP  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if tx_i[i] < 0:
                        file_handle.write("%.9E " % float(tx_i[i]))
                    else:
                        file_handle.write(" %.9E " % float(tx_i[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== TXVAR.EXP section =====
                file_handle.write("\n \n")
                file_handle.write(">TXVAR.EXP  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if tx_var[i] < 0:
                        file_handle.write("%.9E " % float(tx_var[i]))
                    else:
                        file_handle.write(" %.9E " % float(tx_var[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== TYR.EXP section =====
                file_handle.write("\n \n")
                file_handle.write(">TYR.EXP  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if ty_r[i] < 0:
                        file_handle.write("%.9E " % float(ty_r[i]))
                    else:
                        file_handle.write(" %.9E " % float(ty_r[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== TYI.EXP section =====
                file_handle.write("\n \n")
                file_handle.write(">TYI.EXP  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if ty_i[i] < 0:
                        file_handle.write("%.9E " % float(ty_i[i]))
                    else:
                        file_handle.write(" %.9E " % float(ty_i[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                # ===== TYVAR.EXP section =====
                file_handle.write("\n \n")
                file_handle.write(">TYVAR.EXP  //" + str(np.size(freq)))
                file_handle.write("\n")

                for i in range(np.size(freq)):
                    k = i + 1
                    if ty_var[i] < 0:
                        file_handle.write("%.9E " % float(ty_var[i]))
                    else:
                        file_handle.write(" %.9E " % float(ty_var[i]))
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

            if project_setup['processing_mode'] == 'MT + Tipper' or project_setup['processing_mode'] == 'MT Only':
                if procinfo['remotesite'] is not None:
                    # ===== cohEx section =====
                    file_handle.write("\n \n")
                    file_handle.write(
                        ">COH  MEAS1=" + "101.001" + " MEAS2=" + "107.001" + "  ROT=NORTH //" + str(np.size(freq)))
                    file_handle.write("\n")

                    for i in range(np.size(freq)):
                        k = i + 1
                        if coh_ex[i] < 0:
                            file_handle.write("%.9E " % float(coh_ex[i]))
                        else:
                            file_handle.write(" %.9E " % float(coh_ex[i]))
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                    # ===== cohEy section =====
                    file_handle.write("\n \n")
                    file_handle.write(
                        ">COH  MEAS1=" + "102.001" + " MEAS2=" + "106.001" + "  ROT=NORTH //" + str(np.size(freq)))
                    file_handle.write("\n")

                    for i in range(np.size(freq)):
                        k = i + 1
                        if coh_ey[i] < 0:
                            file_handle.write("%.9E " % float(coh_ey[i]))
                        else:
                            file_handle.write(" %.9E " % float(coh_ey[i]))
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")
                else:
                    # ===== cohEx section =====
                    file_handle.write("\n \n")
                    file_handle.write(
                        ">COH  MEAS1=" + "101.001" + " MEAS2=" + "105.001" + "  ROT=NORTH //" + str(np.size(freq)))
                    file_handle.write("\n")

                    for i in range(np.size(freq)):
                        k = i + 1
                        if coh_ex[i] < 0:
                            file_handle.write("%.9E " % float(coh_ex[i]))
                        else:
                            file_handle.write(" %.9E " % float(coh_ex[i]))
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                    # ===== cohEy section =====
                    file_handle.write("\n \n")
                    file_handle.write(
                        ">COH  MEAS1=" + "102.001" + " MEAS2=" + "104.001" + "  ROT=NORTH //" + str(np.size(freq)))
                    file_handle.write("\n")

                    for i in range(np.size(freq)):
                        k = i + 1
                        if coh_ey[i] < 0:
                            file_handle.write("%.9E " % float(coh_ey[i]))
                        else:
                            file_handle.write(" %.9E " % float(coh_ey[i]))
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

            file_handle.write("\n \n")
            file_handle.write(">END")
            file_handle.close()


def read_edi(file_path):
    """
    Docs
    """
    data = {}
    f = open(file_path, 'r')
    edi_cont = f.readlines()

    # # Iterate through each line in the file
    # for line in edi_cont:
    #     line = line.strip()  # Remove leading and trailing whitespaces
    #     if line.startswith('REFLAT'):
    #         reflat = line.split('=')[1]
    #         flat = [float(i) for i in reflat.split(':')]
    #         lat = flat[0] + flat[1] / 60 + flat[2] / (60 * 60)
    #     if line.startswith('REFLONG'):
    #         reflat = line.split('=')[1]
    #         flong = [float(i) for i in reflat.split(':')]
    #         long = flong[0] + flong[1] / 60 + flong[2] / (60 * 60)
    #     if line.startswith('NFREQ'):
    #         nfreq = int(line.split('=')[1])

    # Read Frequencies
    for index, line in enumerate(edi_cont):
        line = line.strip()
        if line.startswith('>FREQ'):
            break
    freqs = []
    for i in range(index + 1, len(edi_cont)):
        line = edi_cont[i]
        try:
            values = [float(val) for val in line.split()]
            freqs.extend(values)
        except ValueError:
            break
    data['freqs'] = freqs

    # Number of coils
    hmeas = 0
    for line in edi_cont:
        if line.startswith('>HMEAS'):
            hmeas += 1

    # =========== MT Impedance =====================
    data['zxx_r'] = read_component('>ZXXR', freqs, edi_cont)
    data['zxx_i'] = read_component('>ZXXI', freqs, edi_cont)
    data['zxx_var'] = read_component('>ZXX.VAR', freqs, edi_cont)
    data['zxy_r'] = read_component('>ZXYR', freqs, edi_cont)
    data['zxy_i'] = read_component('>ZXYI', freqs, edi_cont)
    data['zxy_var'] = read_component('>ZXY.VAR', freqs, edi_cont)
    data['zyx_r'] = read_component('>ZYXR', freqs, edi_cont)
    data['zyx_i'] = read_component('>ZYXI', freqs, edi_cont)
    data['zyx_var'] = read_component('>ZYX.VAR', freqs, edi_cont)
    data['zyy_r'] = read_component('>ZYYR', freqs, edi_cont)
    data['zyy_i'] = read_component('>ZYYI', freqs, edi_cont)
    data['zyy_var'] = read_component('>ZYY.VAR', freqs, edi_cont)

    # =========== Tipper =====================
    data['tzx_r'] = read_component('>TXR.EXP', freqs, edi_cont)
    data['tzx_i'] = read_component('>TXI.EXP', freqs, edi_cont)
    data['tzx_var'] = read_component('>TXVAR.EXP', freqs, edi_cont)
    data['tzy_r'] = read_component('>TYR.EXP', freqs, edi_cont)
    data['tzy_i'] = read_component('>TYI.EXP', freqs, edi_cont)
    data['tzy_var'] = read_component('>TYVAR.EXP', freqs, edi_cont)

    # =========== Coherency =====================
    data['coh_ex'] = read_component('>COH  MEAS1=101.001', freqs, edi_cont)
    data['coh_ey'] = read_component('>COH  MEAS1=102.001', freqs, edi_cont)
    f.close()
    return edi_cont, data, hmeas


def read_component(text, freqs, edi_cont):
    """
    Docs
    """
    data = []
    for index, line in enumerate(edi_cont):
        line = line.strip()
        if line.startswith(text):
            break
    for i in range(index + 1, len(edi_cont)):
        line = edi_cont[i]
        try:
            values = [float(val) for val in line.split()]
            data.extend(values)
        except ValueError:
            break
    if not data:
        data = [math.nan * int(freq) for freq in freqs]
    return data
