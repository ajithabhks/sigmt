import numpy as np
def cohEx(bandavg):
    #Ex predicted
    ExExc = bandavg.get('ExExc')
    ExHxc = bandavg.get('ExHxc')
    ExHyc = bandavg.get('ExHyc')
    HxHxc = bandavg.get('HxHxc')
    HyHxc = bandavg.get('HyHxc')
    HxHyc = bandavg.get('HxHyc')
    HyHyc = bandavg.get('HyHyc')
    Zxx = bandavg.get('Zxx_single')
    Zxy = bandavg.get('Zxy_single')
    ZpZ = Zxx * np.conj(ExHxc) + Zxy * np.conj(ExHyc)
    ZpX = Zxx * HxHxc + Zxy * HyHxc
    ZpY = Zxx * HxHyc + Zxy * HyHyc
    ZpZp = Zxx * np.conj(ZpX) + Zxy * np.conj(ZpY)
    ZpZp = ZpZp * ExExc
    Ccoh = np.empty((ZpZp.shape),dtype=complex)
    for j in range(ZpZp.shape[0]):
        for i in range(ZpZp.shape[1]):
            if abs(ZpZp[j,i])>0:
                Ccoh[j,i] = ZpZ[j,i]/np.sqrt(ZpZp[j,i])
            else:
                Ccoh[j,i] = 1+1j
        cohEx = abs(Ccoh)
        for i in range(cohEx.shape[1]):
            if cohEx[j,i] > 1.0:
                cohEx[j,i] = 1/cohEx[j,i]
    AllcohEx = cohEx
    return AllcohEx
def cohEy(bandavg):
    #Ey predicted
    EyEyc = bandavg.get('EyEyc')
    EyHxc = bandavg.get('EyHxc')
    EyHyc = bandavg.get('EyHyc')
    HxHxc = bandavg.get('HxHxc')
    HyHxc = bandavg.get('HyHxc')
    HxHyc = bandavg.get('HxHyc')
    HyHyc = bandavg.get('HyHyc')
    Zyy = bandavg.get('Zyy_single')
    Zyx = bandavg.get('Zyx_single')
    ZpZ = Zyx * np.conj(EyHxc) + Zyy * np.conj(EyHyc)
    ZpX = Zyx * HxHxc + Zyy * HyHxc
    ZpY = Zyx * HxHyc + Zyy * HyHyc
    ZpZp = Zyx * np.conj(ZpX) + Zyy * np.conj(ZpY)
    ZpZp = ZpZp * EyEyc
    Ccoh = np.empty((ZpZp.shape),dtype=complex)
    for j in range(ZpZp.shape[0]):
        for i in range(ZpZp.shape[1]):
            if abs(ZpZp[j,i])>0:
                Ccoh[j,i] = ZpZ[j,i]/np.sqrt(ZpZp[j,i])
            else:
                Ccoh[j,i] = 1+1j
        cohEy = abs(Ccoh)
        #cohEy = 1-cohEy
        for i in range(cohEy.shape[1]):
            if cohEy[j,i] > 1.0:
                cohEy[j,i] = 1/cohEy[j,i]
    AllcohEy = cohEy
    return AllcohEy
