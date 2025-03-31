"""Optimises stage sizes for maximum dV

Rockets are represented by n-stages
- auxilliary mass
- propulsion system (with isp, wet and dry mass) <- this is what is optimised for, auxilliary mass cannot be changed
- max burn times, so we don't overstress an engine

To scan efficiently over the available parameter space, we treat this as a sphere traversal problem. We don't guess the stage masses directly, but 
a position on an n-sphere of given radius using the 'surface' degrees of freedom. This guess corresponds to a rocket of some total mass with stage masses
that add to the correct value.
"""

import numpy as np
from math import cos, pi, sqrt, log
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# TODO: estimate a function of propLambda vs propulsion system mass, given the propulsion system type: (solid, liquid, hybrid)
# TODO: account for tank residuals, some percentage in each tank

def angs2masses(totalMass:float, angs:list[float]) -> list[float]:
    """
    Converts spherical coordinates to rocket stage masses that sum to a specified total mass.

    The angles are angles used to specify the position, p, on an n-sphere's surface. The rocket's mass is 
    fed in as well, and converted to a representative radius, r.

    As the stage's masses must sum to some total, mTotal, we can use the fact that on an n-sphere:
    sqrt(sum(coords**2)) = r.

    Since we want some value for each stage, we take each stage as coords[i]**2, and r**2 is the total vehicle mass.
    """
    
    r = sqrt(totalMass)
    angs = np.array(angs, float)

    coords = np.zeros((angs.size + 1), float) # cartesian coordinates

    for i in range(0, angs.size):
        coords[i] = r * np.prod(np.sin(angs[0:i])) * cos(angs[i])

    coords[-1] = r * np.prod(np.sin(angs))

    # convert to stage masses:
    return coords**2



class Stage():
    """
    Represents one stage of a launch vehicle. Each stage consists of a propulsion system and some extra mass.
    The propulsion systems is made up from inert mass and propellant. The propellant is fully used during the stage's lifetime.
    The stage also has an associated thrust, and in future, maximum burn time (which couples with thrust or maximum propellant mass)

    Planned features:
    - automatic estimate of structural efficiency using scaling laws and a named reference stages
  
    *Constraints*
    Constraints are optional, and can be used to ensure that the 'optimum' vehicle found takes other factors into consideration. For example, perhaps the engine on stage 2 is only rated to run for 60 seconds.

    Available stage constraints:
    - mass limits: [min, max]
    - burn time [min, max]
    - acceleration [min, max]

    if your constraint does not have a limit, then put either -float('inf') (or 0 if more appropriate) or float('inf') for min or max respectively - i.e. constraints['burn time'] = [0, float('inf')]
    if the constraint is a single number, make both the minimum and maximum the same value (or make it a single number and I get the type?)

    Reference Stages and estimated structural efficiency:
    - uses scaling laws and reference stages to predict the structural efficiency of a given stage.
    """

    def __init__(self, totalMass:float, isp:float, massFlow:float, propLambda:float, extraMass:float, canVary:bool, constraints:dict={}):

        self.totalMass = totalMass # total mass of the stage
        self.isp = isp
        self.massFlow = massFlow
        self.propLambda = propLambda
        self.extraMass = extraMass
        self.canVary = canVary # if false, do not change this stage when optimising
        self.constraints = constraints

        self.calculate()

    
    def calculate(self, newTotalMass:float=None, newIsp:float=None, newMassFlow:float=None, newPropLambda:float=None, newExtraMass:float=None):
        """Get useful information about the stage given inputs"""

        if newTotalMass is not None:
            self.totalMass = newTotalMass

        if newIsp is not None:
            self.isp = newIsp

        if newMassFlow is not None:
            self.massFlow = newMassFlow

        if newPropLambda is not None:
            self.propLambda = newPropLambda

        if newExtraMass is not None:
            self.extraMass = newExtraMass

        self.mInit = self.totalMass
        self.mPropellant = (self.totalMass - self.extraMass) * self.propLambda
        self.mFinal = self.mInit - self.mPropellant
        self.burnTime = self.mPropellant / self.massFlow
        self.thrust = self.isp * self.massFlow * 9.81



class Rocket():
    """
    Rockets consist of n-stages as well as a payload.
    """

    def __init__(self, stages:list[Stage], payload:float):

        self.stages = stages
        self.payload = payload
        self.canVary, self.freeStages = self.getVariableStages()


    def deltaV(self):
        # get the delta v contribution by each stage

        deltaVs = np.zeros((len(self.stages)), float)
        initMasses = np.zeros((len(self.stages)), float)
    
        for i in range(0, len(self.stages)):
            initMasses[i] = self.stages[i].mInit

        for i in range(0, len(self.stages)):

            m0 = np.sum(initMasses[i:]) + self.payload
            mf = m0 - self.stages[i].mPropellant
            deltaVs[i] = self.stages[i].isp * 9.81 * log(m0/mf)

        return deltaVs

    
    def calcAccelerations(self):
        """Return a list of maximum accelerations by stage"""
        
        maxAccelerations = []

        for i in range(0, len(self.stages)):
            # we assume a constant acceleration:
            mFinal = -self.stages[i].mPropellant

            for upperStage in self.stages[i:]:
                mFinal += upperStage.mInit

            maxAccelerations.append(self.stages[i].thrust / mFinal)

        return maxAccelerations

    
    def getVariableStages(self):

        variableStages = []
        freeStages = []

        for i in range(0, len(self.stages)):
            variableStages.append(self.stages[i].canVary)
            
            if self.stages[i].canVary:
                freeStages.append(i)

        return variableStages, freeStages 

    
    def changeStage(self, stageNum:int, newTotalMass:float=None, newIsp:float=None, newPropLambda:float=None, newExtraMass:float=None):
 
        self.stages[stageNum].calculate(newTotalMass, newIsp, newPropLambda, newExtraMass)

    
    def addStage(self, newStage:Stage, num:int):
        """Inserts a new stage at the desired position (1-indexed)"""
        self.stages.insert(num, newStage)

    
    def removeStage(self, number:int):
        """Removes the stage with the corresponding number (1-indexed)"""
        self.stages.pop(number - 1)


    def scoreStages(self) -> bool:
        """
        Checks every stage to see if their constraints were satisfied.
        
        Each stage is only checked on their specific constraints, returning a dict of bools to store in a dict or list.
        """

        self.stageScores = []
        totalDeviation = 0

        for i in range(0, len(self.stages)):
            
            stage = self.stages[i]
            score = {}

            massLimits = stage.constraints.get('mass')
            accLimits = stage.constraints.get('acceleration')
            timeLimits = stage.constraints.get('time')

            if massLimits is not None: # the number of kilograms the stage is outside of bounds.
                
                score['mass'] = max(max(0, massLimits[0] - stage.mInit), max(0, stage.mInit - massLimits[1]))
                totalDeviation += score['mass']

            if accLimits is not None: # the amount of m/s^2 the acceleration lies outside of bounds
                
                massInit = 0
                for upperStage in self.stages[i:]:
                    massInit += upperStage.mInit

                massFinal = massInit - stage.mPropellant

                minAcc = stage.thrust / massInit
                maxAcc = stage.thrust / massFinal

                score['acceleration'] = max(max(0, accLimits[0] - minAcc), max(0, maxAcc - accLimits[1]))
                totalDeviation += score['acceleration']

            if timeLimits is not None: # seconds by which the stage burn time lies outside of bounds
                
                #print(f"Propellant mass: {stage.mPropellant}, Mass flow-rate: {stage.massFlow}, Burn time: {stage.burnTime}")
                score['time'] = max(max(0, timeLimits[0] - stage.burnTime), max(0, stage.burnTime - timeLimits[1]))
                totalDeviation += score['time']

            self.stageScores.append(score)

        # did any of the checks fail?
        return totalDeviation == 0.0



def angs2dv(angs:list[float], rocketMass:float, rocket:Rocket) -> float:
    """This function takes in spherical coordinates and outputs the rocket's total delta-v (* -1 for use in scipy minimise functions)"""

    # calculate the available mass to allocate:
    fixedMass = 0
    for stage in rocket.stages:
        fixedMass += stage.totalMass if not stage.canVary else 0

    freeMass = rocketMass - rocket.payload - fixedMass

    stageMasses = angs2masses(totalMass=freeMass, angs=angs)

    for i in range(0, len(rocket.freeStages)):
        rocket.changeStage(stageNum=rocket.freeStages[i], newTotalMass=stageMasses[i])

    rocket.scoreStages()

    massMultiplier = 1000
    accMultiplier = 1000
    timeMultiplier = 1000

    # we may be able to handle this better but for now, just add up all deviations by type:
    massPenalty = 0
    accPenalty = 0
    timePenalty = 0

    for score in rocket.stageScores:

        massScore = score.get('mass')
        accScore = score.get('acceleration')
        timeScore = score.get('time')

        massPenalty += 0 if massScore is None else massMultiplier * massScore
        accPenalty += 0 if accScore is None else accMultiplier * accScore
        timePenalty += 0 if timeScore is None else timeMultiplier * timeScore

    return -np.sum(rocket.deltaV()) + massPenalty + accPenalty + timePenalty



def optimise(rocket:Rocket, rocketMass:float) -> Rocket:
    """Calculates a rocket with a given total mass that delivers the highest possible delta-v according to design constraints."""

    # set input limits:
    angLims = [-pi/2, pi/2]
    bounds = []

    for i in range(0, len(rocket.freeStages)-1):#
        bounds.append(angLims)

    bounds = tuple(bounds)
    
    x0 = (0.5 * np.ones((len(rocket.freeStages) - 1), float)).tolist()
  
    minimize(fun=angs2dv, x0=x0, args=(rocketMass, rocket), method='Nelder-Mead', bounds=bounds, tol=1e-9)

    return rocket



def sweepMass(rocket:Rocket, limits=list[float]):
    """Visualise how vehicle performance changes by varying different properties of the launch vehicle"""

    # let's sweep the vehicle total mass:
    rocketMass = np.linspace(limits[0], limits[1], 50)
    currentMass = 10000 # use this to see how far from optimum a certain design is, given its current mass

    # Moonlight ICBM: when playing KSP, propLambda is just the tank utilisation, everything else is either extra mass or payload
    # stage1 = Stage(totalMass=0, isp=215, massFlow=53.9, propLambda=0.918, extraMass=574, canVary=True, constraints={'time':[0,103]}) # MOONLIGHT BOOSTER
    # stage2 = Stage(totalMass=0, isp=236, massFlow=8.64, propLambda=0.785, extraMass=76, canVary=True, constraints={'time':[0,60]}) # MOONLIGHT 2nd Stage
    # stage3 = Stage(totalMass=0, isp=228, massFlow=6.15, propLambda=0.764, extraMass=25, canVary=True, constraints={'time':[0,50]}) # MOTH M


    dv = np.zeros((rocketMass.size, len(rocket.stages)), float)
    maxAccelerations = np.zeros((rocketMass.size, len(rocket.stages)), float)
    stageMasses = np.zeros((rocketMass.size, len(rocket.stages)), float)
    stageTimes = np.zeros((rocketMass.size, len(rocket.stages)), float)
    propMasses = np.zeros((rocketMass.size, len(rocket.stages)), float)

    for i in range(0, rocketMass.size):
        
        optimRocket = optimise(rocket, rocketMass[i])

        for j in range(0, len(optimRocket.stages)):
            stageMasses[i,j] = optimRocket.stages[j].mInit
            stageTimes[i,j] = optimRocket.stages[j].burnTime
            propMasses[i,j] = optimRocket.stages[j].mPropellant

        dv[i,:] = optimRocket.deltaV()
        maxAccelerations[i,:] = optimRocket.calcAccelerations()

    currentMassX = [currentMass, currentMass]
    currentMassY = [np.max(np.sum(dv, axis=1)), np.min(np.sum(dv, axis=1))]

    fig, axs = plt.subplots(2, 3)

    axs[0,0].plot(rocketMass, np.sum(dv, axis=1), '-')
    axs[0,0].plot(currentMassX, currentMassY, '--r')
    axs[0,0].set_xlabel('rocket mass, kg')
    axs[0,0].set_ylabel('optimum Δv, m/s')
    axs[0,0].grid()

    stageKey = []

    for i in range(0, len(rocket.stages)):
        axs[1,0].plot(rocketMass, dv[:,i])
        axs[0,1].plot(rocketMass, stageMasses[:,i])
        axs[1,1].plot(rocketMass, maxAccelerations[:,i])
        axs[0,2].plot(rocketMass, stageTimes[:,i])
        axs[1,2].plot(rocketMass, propMasses[:,i])
        stageKey.append('stage' + str(i + 1))

    axs[1,0].set_xlabel('rocket mass, kg')
    axs[1,0].set_ylabel('Δv, m/s')
    axs[1,0].legend(stageKey)
    axs[1,0].grid()

    axs[0,1].set_xlabel('rocket mass, kg')
    axs[0,1].set_ylabel('stage mass, kg')
    axs[0,1].legend(stageKey)
    axs[0,1].grid()
    
    axs[1,1].set_xlabel('rocket mass, kg')
    axs[1,1].set_ylabel('max acceleration, m/s^2')
    axs[1,1].legend(stageKey)
    axs[1,1].grid()

    axs[0,2].set_xlabel('rocket mass, kg')
    axs[0,2].set_ylabel('stage burn time, s')
    axs[0,2].legend(stageKey)
    axs[0,2].grid()

    axs[1,2].set_xlabel('rocket mass, kg')
    axs[1,2].set_ylabel('propellant mass, kg')
    axs[1,2].legend(stageKey)
    axs[1,2].grid()

    plt.show()



def targetDv(rocket:Rocket, dvTarget:float, massLims:list):
    """Calculate the lightest rocket that could achieve the target ideal delta-v"""

    # method 1: sweep over mass range, get crossing points of function, then solve using minimize or shooting method over this small range
    dm = 50
    rocketMass = np.arange(massLims[0], massLims[1], dm)
    dv = [0]
    i = 0

    while(dv[-1] <= dvTarget and i < rocketMass.size):
        
        finalRocket = optimise(rocket, rocketMass[i])
        dv.append(sum(finalRocket.deltaV()))
        i += 1

    if i == rocketMass.size:
        print("No vehicle within mass range achieved target dv, please adjust settings and try again.")
        return
      
    else:
        # there is a crossing point: solve for the mass between rocketMass[i-1] and rocketMass[i-2]
        mass = rocketMass[i-2:i].tolist()

        tol = 1e-3
        err = [2*tol, 2*tol]
        kRelax = 1

        while(abs(err[1]) > tol):
            # Newton's method for rocket mass
            err[0] = sum(optimise(rocket, mass[0]).deltaV()) - dvTarget
            err[1] = sum(optimise(rocket, mass[1]).deltaV()) - dvTarget

            newMass = mass[0] + (err[0] / (err[0] - err[1])) * (mass[1] - mass[0])

            mass[0] = mass[1]
            mass[1] = mass[1] * (1-kRelax) + newMass * kRelax # NOTE: we can apply a relaxation factor if this is unstable

        # get the outputs for the rocket
        finalRocket = optimise(rocket, mass[1])
        finalDv = sum(finalRocket.deltaV())
        isValid = finalRocket.scoreStages()

        print(f"final rocket mass: {mass[1]} kg, delta-v: {finalDv} m/s (error = {err[1]})")
        print(f"checks passed: {isValid}")

        # output the stage breakdown for this rocket:
        for j in range(0, len(rocket.stages)):
            print(f"Stage {j}:\n=========================================================")
            print(f"Total Mass:\t\t{'%.3f' % rocket.stages[j].totalMass} kg")
            print(f"Propellant Mass:\t{'%.3f' % rocket.stages[j].mPropellant} kg")
            print(f"Final Mass:\t\t{'%.3f' % rocket.stages[j].mFinal} kg")
            print(f"Required burn time:\t{'%.3f' % rocket.stages[j].burnTime} s")
            print(f"delta-v:\t\t{'%.3f' % rocket.deltaV()[j]} m/s")
            print('\n')

    fig, ax = plt.subplots()
    ax.plot(rocketMass[:i], dv[:i], '-k')
    ax.set_xlabel('rocket mass, kg')
    ax.set_ylabel('total dV, m/s')

    plt.show()



def main():

    # build the rocket:
    stage1 = Stage(totalMass=4332, isp=290.65, massFlow=50, propLambda=0.917, extraMass=40, canVary=False, constraints={})
    stage2 = Stage(totalMass=0, isp=340, massFlow=24.662, propLambda=0.9, extraMass=40, canVary=True, constraints={})
    stage3 = Stage(totalMass=0, isp=340, massFlow=4.932, propLambda=0.9, extraMass=20, canVary=True, constraints={})

    stages = [stage1, stage2, stage3]
    rocket = Rocket(stages=stages, payload=250)

    # sweepMass(rocket=rocket, limits=[6000, 10000])
    targetDv(rocket=rocket, dvTarget=8500, massLims=[5000, 50000])
    #optimise(rocket=rocket, rocketMass=8500)
    

main()