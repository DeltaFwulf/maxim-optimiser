import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from math import copysign

from rocket import Rocket, Stage

def main():

    planet = {}
    planet.update({'r':6378e3})
    planet.update({'mu':3.986e14})

    payload = 100.0
    
    tf = 300.0
    dt = 0.1
    t = np.arange(0, tf, dt)

    stages = []
    stages.append(Stage(mTotal=1000, propFrac=0.8, isp=250, mDot=30, canVary=False, cdA=0.2))
    stages.append(Stage(mTotal=500, propFrac=0.7, isp=300, mDot=2, canVary=False, cdA=0.1))
    rocketA = Rocket(stages)
    rocketA.tVals = [1.0, 1.0]
    rocketA.tTimes = [0.0, tf]
    rocketA.name = 'full throttle'

    stages = []
    stages.append(Stage(mTotal=1000, propFrac=0.8, isp=250, mDot=30, canVary=False, cdA=0.2))
    stages.append(Stage(mTotal=500, propFrac=0.7, isp=300, mDot=2, canVary=False, cdA=0.1))
    rocketB = Rocket(stages)
    rocketB.tVals = [0.8, 0.8]
    rocketB.tTimes = [0.0, tf]
    rocketB.name = '80 percent throttle'

    stages = []
    stages.append(Stage(mTotal=1000, propFrac=0.8, isp=250, mDot=30, canVary=False, cdA=0.2))
    stages.append(Stage(mTotal=500, propFrac=0.7, isp=300, mDot=2, canVary=False, cdA=0.1))
    rocketC = Rocket(stages)
    rocketC.tVals = [1.0, 1.0, 0.5, 0.5, 1.0, 1.0]
    rocketC.tTimes = [0.0, 12.0, 13.0, 21.0, 23.0, tf]
    rocketC.name = 'short damping'


    rockets = [rocketA, rocketB, rocketC]

    z0 = 0.0
    v0 = 0.0
    
    fig, ax = plt.subplots(3, 1)

    ax[0].set_xlabel('time, s')
    ax[0].set_ylabel('z, m')

    ax[1].set_xlabel('time, s')
    ax[1].set_ylabel('v, m/s')
    
    ax[2].set_xlabel('time, s')
    ax[2].set_ylabel('q, Pa')

    for roc in rockets:

        z, v, q = maxQ(z0, v0, tf, dt, roc, planet, payload)

        ax[0].plot(t, z, '-', label=roc.name)
        ax[1].plot(t, v, '-', label=roc.name)
        ax[2].plot(t, q, '-', label=roc.name)
        
    ax[0].legend()
    #ax[1].legend(loc='lower right')
    #ax[2].legend()
    plt.show()



def maxQ(z0, v0, tf, dt, rocket, planet, payload):

    t = np.arange(0, tf, dt)
    z = np.zeros_like(t)
    v = np.zeros_like(t)
    q = np.zeros_like(t)

    # Initialise arrays:
    z[0] = z0
    v[0] = v0

    tVals = rocket.tVals
    tTimes = rocket.tTimes

    # run the simulation until time is up 
    for i in range(1, t.size):

        throttle = np.interp(t[i-1], tTimes, tVals)

        params = {}
        params.update({'planet':planet})
        params.update({'payload':payload})
        params.update({'rocket':deepcopy(rocket)})
        params.update({'throttle':throttle})

        # rk4 update state
        xIn = np.array([z[i-1], v[i-1]], float)
        k1 = diff(xIn, params)

        # update the vehicle state (assuming constant thrust over the half timestep)
        rocket.simulate(dt / 2, isBurning=True, throttle=throttle)

        if rocket.getCurrentStage().mp == 0:
            rocket.advanceStage()

        throttle = np.interp(t[i-1] + dt / 2, tTimes, tVals)

        params.update({'rocket':deepcopy(rocket)})
        params.update({'throttle':throttle})

        k2 = diff(xIn + k1*dt / 2, params)
        k3 = diff(xIn + k2*dt / 2, params)

        rocket.simulate(dt / 2, isBurning=True, throttle=throttle) # gives the final state of the vehicle

        if rocket.getCurrentStage().mp == 0:
            rocket.advanceStage()

        throttle = np.interp(t[i-1] + dt, tTimes, tVals)

        params.update({'rocket':deepcopy(rocket)})
        params.update({'throttle':throttle})

        k4 = diff(xIn + k3*dt, params)
        xOut = xIn + (dt / 6)*(k1 + 2*(k2 + k3) + k4)

        z[i] = xOut[0] if xOut[0] > 0.0 else 0.0
        v[i] = xOut[1] if xOut[0] > 0.0 else 0.0    
        q[i] = 0.5*coesa76(z[i])['density']*v[i]**2

    return z, v, q

     

def diff(x:np.array, params:dict) -> np.array:

    roc = params['rocket']
    throttle = params['throttle']
    planet = params['planet']
    payload = params['payload']
    
    T = roc.getThrust(throttle)
    rho = coesa76(x[0])['density']
    D = -copysign(1.0, x[1]) * 0.5*rho*x[1]**2 * roc.getCurrentStage().cdA # opposes direction of motion
    m = roc.getMass() + payload
    g = planet['mu'] / (planet['r'] + x[0])**2
    
    a = (T + D - m*g) / m

    return np.array([x[1], a], float)


def coesa76(z:float) -> dict:

    """interpolates COESA76 data for temperature, density, and pressure"""

    zArr = [0, 11019, 20063, 32162, 47350, 51413, 71802, 86000, 91000, 110000, 120000, 500000, 1000000]
    Tarr = [288.15, 216.65, 216.65, 288.65, 270.65, 270.65, 214.65, 186.87, 186.87, 240.00, 360.00, 999.24, 1000.00]
    pArr = [101325, 22632, 5474.8, 868.01, 110.90, 66.938, 3.9564, 0.37338, 0.15381, 0.0071042, 0.0025382, 0.00000030236, 0.0075138]
    rhoArr = [1.225, 3.6392e-1, 8.8035e-2, 1.3225e-2, 1.4275e-3, 8.6160e-4, 6.4211e-5, 6.958e-6, 2.860e-6, 9.708e-8, 2.222e-8, 5.215e-13, 3.561e-15]
    
    T = np.interp(z, zArr, Tarr)
    p = np.interp(z, zArr, pArr)
    rho = np.interp(z, zArr, rhoArr)

    atmo = {}
    atmo.update({'temperature':T})
    atmo.update({'pressure':p})
    atmo.update({'density':rho})

    return atmo


main()