import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from math import exp

from rocket import Stage, Rocket
from optimiser import maxPayload



def maxPayloadVsMass():
    """given a rocket and stages with associated scaling constants, estimate maximum payload vs rocket total mass for a given dv requirement""" 
    
    dv = 8750.0
    TWR = 2.0

    stages = []
    stages.append(Stage(mTotal=4292+60, propFrac=0.915, isp=296, mDot=55, canVary=False, fixedMass=60, constraints={}))
    stages.append(Stage(mTotal=873+50, propFrac=0.876, isp=287, mDot=11, canVary=False, fixedMass=50, constraints={'acceleration':[20, 60]}))
    stages.append(Stage(mTotal=925, propFrac=0.79, isp=373, mDot=4.1, canVary=True, fixedMass=50, constraints={}))
    rocket = Rocket(stages)

    propellants = []
    propellants.append({'rho':942.424, 'isp':372, 'name':'Propane', 'color':'red'})
    propellants.append({'rho':1011.01, 'isp':366, 'name':'RP-1', 'color':'black'})
    propellants.append({'rho':830, 'isp':378, 'name':'Methane', 'color':'blue'})

    mp = np.arange(50, 4000, 10)
    mpl = np.zeros_like(mp, float)
    dvs = np.zeros((3, mp.size), float)
    p = np.zeros_like(mp, float)
    pMin = np.zeros_like(mp, float)
    stageMasses = np.zeros_like(dvs, float)
    mtot = np.zeros_like(mp, float)

    fiVpRef = 0.116
    cf = 0.0101
    ampRef = 0.0839
    mpRef = 638.8
    rhoRef = 942.424
    
    fig0, ax0 = plt.subplots()
    fig1, ax1 = plt.subplots()

    bigMax = 0.0

    for prop in propellants:

        rocket.stages[-1].update({'isp':prop['isp']})
        cmpRef = ampRef * (prop['rho'] / rhoRef)**-0.75
        firstPass = -1

        for i in range(mp.size):

            # calculate stage 3 total mass
            fi = fiVpRef*(rhoRef / prop['rho']) + cf*TWR + cmpRef*(mp[i] / mpRef)**(-2/3)
            pProp = 1 - fi
            mStage = mp[i] / pProp + rocket.stages[-1].fixedMass
            rocket.stages[-1].update({'total_mass':mStage, 'prop_fraction':pProp})
            p[i] = rocket.stages[-1].mp / rocket.stages[-1].mtot
            
            # vehicle performance
            rocketOut, mpl[i] = maxPayload(rocket, rocket.getMass(), dv)
            mtot[i] = rocket.getMass()
            stageDvs = rocketOut.getDv(mpl[i])
            pMin[i] = 1 - exp(-stageDvs[-1] / (rocketOut.stages[-1].isp * 9.80665))
            
            for s in range(len(rocket.stages)):
                dvs[s, i] = stageDvs[s]
                stageMasses[s, i] = rocketOut.stages[s].mtot

            if mpl[i] > 0 and firstPass < 0:
                firstPass = i

            elif mpl[i] <= 0 and firstPass >= 0:
                break

        if firstPass < 0:
            print(f"{prop['name']} did not pass for any propellant masses")
            continue

        mpl_max = np.max(mpl[firstPass:i])
        i_max = np.where(mpl == mpl_max)[0]
        bigMax = mpl_max if mpl_max > bigMax else bigMax

        ax0.plot(mtot[firstPass:i], mpl[firstPass:i], '-', color=prop['color'], label=prop['name'])
        ax0.plot(mtot[i_max], mpl_max, 'o', color=prop['color'], label='_nolegend_')
        ax0.plot([mtot[i_max], mtot[i_max]], [0, mpl_max], '--', linewidth=0.8, color=prop['color'], label='_nolegend_')

        ax1.plot(mtot[firstPass:i], p[firstPass:i], '-', linewidth=1, color=prop['color'], label=f"expected {prop['name']}")
        ax1.plot(mtot[firstPass:i], pMin[firstPass:i], '--', linewidth=0.8, color=prop['color'], label=f"min {prop['name']}")
        
    ax0.set_title(f"Maximum Payload vs Launch Mass for Δv = {'%i' % dv} m/s")
    ax0.set_xlabel("Rocket mass, kg")
    ax0.set_ylabel("Max. payload, kg")
    ax0.xaxis.set_major_locator(ticker.MultipleLocator(500))
    ax0.xaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax0.set_ylim([0, bigMax*1.05])
    ax0.legend()
    ax0.grid(True)

    ax1.set_title(f"Final stage λ vs Launch Mass for Δv = {'%i' % dv} m/s")
    ax1.set_xlabel("Launch mass, kg")
    ax1.set_ylabel("λ")
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(500))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax1.legend()
    ax1.grid(True)

    plt.show()



maxPayloadVsMass()