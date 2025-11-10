import numpy as np
from math import log



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

    def __init__(self, stages:list[Stage]):

        self.stages = stages
        self.canVary, self.freeStages = self.getVariableStages()


    def getMass(self):

        m = 0

        for stage in self.stages:
            m += stage.totalMass

        return m
    

    def deltaV(self, payload:float) -> list['float']:
        # get the delta v contribution by each stage

        stageDv = np.zeros((len(self.stages)), float)
        initMasses = np.zeros((len(self.stages)), float)
    
        for i in range(0, len(self.stages)):
            initMasses[i] = self.stages[i].mInit

        for i in range(0, len(self.stages)):

            m0 = np.sum(initMasses[i:]) + payload
            mf = m0 - self.stages[i].mPropellant
            stageDv[i] = self.stages[i].isp * 9.81 * log(m0/mf)

        return stageDv

    
    def calcAccelerations(self, payload:float) -> list['float']:
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