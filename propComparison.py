import numpy as np
import matplotlib.pyplot as plt
from math import log


def estimateFraction(techLevel:str, mp:float, rho:float, TWR:float):

    # fiVpRef = {'low':0.03, 'med':0.022, 'high':0.015}
    # cf = {'low':0.04, 'med':0.031, 'high':0.025}
    # ampRef = {'low':0.186, 'med':0.127, 'high':0.077}

    fiVpRef = {'med':0.116}
    cf = {'med':0.0101}
    ampRef = {'med':0.0839}

    #mpRef = 4536
    #rhoRef = 1011.022

    mpRef = 868
    rhoRef = 942.424

    cmpRef = ampRef[techLevel] * (rho / rhoRef)**-0.75

    fi = fiVpRef[techLevel]*(rhoRef / rho) + cf[techLevel]*TWR + cmpRef*(mp / mpRef)**(-2/3)
    return 1 - fi



def comparePropellants():

    rp1 = {'rho':1011.022, 'isp':366, 'name':'RP1/LO2'}
    lh2 = {'rho':325.238, 'isp':456, 'name':'LH2/LO2'}
    lprop = {'rho':942.424, 'isp':372, 'name':'Propane/LO2'}
    ch4 = {'rho':830.0, 'isp':378, 'name':'Methane/LO2'}

    props = [rp1, ch4, lprop]

    mp = np.linspace(4000, 6000, 20)
    mpl = 2000.0

    techLevel = 'med'
    TWR = 2.0

    fig, axP = plt.subplots()
    fig1, axM = plt.subplots()
    fig2, axL = plt.subplots()

    dvMin = 0.0
    dvMax = 0.0

    for prop in props:

        dv = np.zeros_like(mp)
        pFrac = np.zeros_like(mp)

        for i in range(mp.size):
            pFrac[i] = estimateFraction(techLevel, mp[i], prop['rho'], TWR)
 
            mi = mp[i]*(1 - pFrac[i]) / pFrac[i]
            m0 = mp[i] / pFrac[i] + mpl
            mf = mi + mpl

            dv[i] = prop['isp'] * 9.80665 * log(m0 / mf)

        axP.plot(mp, dv, '-', label=prop['name'])
        axM.plot(mp/pFrac, dv, '-', label=prop['name'])
        axL.plot(mp, pFrac, '-', label=prop['name'])

        if np.min(dv) < dvMin:
            dvMin = np.min(dv)

        if np.max(dv) > dvMax:
            dvMax = np.max(dv)

    axP.set_xlabel('propellant mass, kg')
    axP.set_ylabel('delta v, m/s')
    axP.set_title(f'delta v vs propellant mass, payload mass = {mpl}, TWR = {TWR}')
    axP.legend()
    axP.grid(True)
    
    axM.set_xlabel('stage wet mass, kg')
    axM.set_ylabel('delta v, m/s')
    axM.set_title(f'delta v vs stage wet mass, payload mass = {mpl}, TWR = {TWR}')
    axM.legend()
    axM.grid(True)

    axL.set_xlabel("propellant mass, kg")
    axL.set_ylabel("propellant mass fraction")
    axL.set_title(f"propellant mass fraction vs propellant mass, TWR = {TWR}")
    axL.legend()
    axL.grid(True)

    plt.show()


comparePropellants()