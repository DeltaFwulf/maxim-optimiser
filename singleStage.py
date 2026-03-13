"""Plot performance trends for single-stage rockets"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp, log


def minFracVsDv(minDv:float, maxDv:float, isps:list['float']) -> None:
    """given a set of isps, draw plots showing how lambda must increase to meet the dv requirements"""
    
    dv = np.linspace(minDv, maxDv, 20)

    fig, ax = plt.subplots()

    for i in range(len(isps)):

        fracs = np.zeros_like(dv)
        for v in range(dv.size):
            fracs[v] = 1 - exp(-dv[v] / (isps[i] * 9.80665))

        ax.plot(dv, fracs, '-', label=str(isps[i]) + ' s')

    ax.legend()
    ax.set_xlabel("dv, m/s")
    ax.set_ylabel("min propellant fraction")
    ax.set_title("Single stage, minimum propellant fraction vs. dv for different isps")
    ax.grid(True)

    plt.show()


def minIspVsLambda(minDv:float, maxDv:float, fracs:list['float']) -> None:
    """Over a range of dv requirements, generate curves for minimum isp given fixed propellant fractions"""
    
    dv = np.linspace(minDv, maxDv, 20)

    fig, ax = plt.subplots()

    for l in range(len(fracs)):

        isps = np.zeros_like(dv)
        for v in range(dv.size):
            isps[v] = dv[v] / (9.80665 * log(1 / (1 - fracs[l])))

        ax.plot(dv, isps, '-', label=str(fracs[l]))

    ax.legend()
    ax.set_xlabel("dv, m/s")
    ax.set_ylabel("min isp, s")
    ax.set_title("Single stage, minimum isp vs dv for fixed propellant fractions")
    ax.grid(True)

    plt.show()



def main():

    minDv = 1000.0
    maxDv = 8000.0

    isps = [270.0, 300.0, 330.0, 360.0, 390.0]
    fracs = [0.7, 0.75, 0.8, 0.85, 0.9]
    #minFracVsDv(minDv, maxDv, isps)

    minIspVsLambda(minDv, maxDv, fracs)

main()