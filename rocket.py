import numpy as np
from math import log, copysign



class Stage():

    # A stage is some object with propellant, inert structure, and a propulsion system

    def __init__(self, mTotal:float, propFrac:float, isp:float, mDot:float, canVary:bool, fixedMass:float=0.0, constraints:dict={}, cdA:float=0.0) -> None:

        self.canVary = canVary
        self.constraints = constraints
        self.cdA = cdA

        params = {}
        params.update({'total_mass':mTotal})
        params.update({'prop_fraction':propFrac})
        params.update({'isp':isp})
        params.update({'mass_flow':mDot})
        params.update({'fixed_mass':fixedMass})
        self.update(params)

    def update(self, params:dict) -> None:

        if params.get('total_mass') is not None:
            self.mtot = params['total_mass']

        if params.get('prop_fraction') is not None:
            self.pFrac = params['prop_fraction']

        if params.get('isp') is not None:
            self.isp = params['isp']

        if params.get('mass_flow') is not None:
            self.mDot = params['mass_flow']

        if params.get('fixed_mass') is not None:
            self.fixedMass = params['fixed_mass']

        # account for fixed mass
        freeMass = self.mtot - self.fixedMass
        
        self.mp = freeMass * self.pFrac
        self.mi = self.mtot - self.mp
        self.thrust = self.isp * 9.80665 * self.mDot
        self.burnTime = self.mp / self.mDot

    def simulate(self, dt, isBurning, throttle) -> None:
        dm = dt * (self.mDot if isBurning else 0) * throttle
        self.mp = max(0.0, self.mp - dm)
        self.mtot = self.mi + self.mp


class ScaledStage():

    def __init__(self, ref:dict, mTot:float, isp:float, mdot:float, rhoProp:float, canVary:bool):
       
        self.fiv_ref = ref['fi_v']
        self.rho_ref = ref['rho']
        self.cf = ref['cf']
        self.cmp = ref['cmp']
        self.mp_ref = ref['mp']

        self.mtot = mTot
        self.isp = isp
        self.mDot = mdot
        self.rho = rhoProp
        self.thrust = mdot*isp*9.80665

        self.canVary = canVary

        # calculate the stage at the set total mass
        self.solveTotalMass()


    def calcMass(self, mp:float):
        """Calculates stage properties given a propellant mass"""

        fi_v = self.fiv_ref*(self.rho_ref / self.rho)
        fi_f = self.cf*self.mDot*self.isp / mp
        fi_m = self.cmp*(self.mp_ref / mp)**(2/3)

        fi = fi_v + fi_f + fi_m
        m = mp / (1-fi)

        return m

    def solveTotalMass(self):
        """Iterates for a given total stage mass"""

        m = self.mtot
        mp = 0
        dmp = 10
        mOut = -1

        # low values lead to negative outputs. Get to the positive region
        while(mOut < 0):
            mp += dmp
            mOut = self.calcMass(mp)

        err0 = mOut - m
        err = err0

        # Now, get the first crossing point
        while copysign(1.0, err) == copysign(1.0, err0):
            mpLast = mp
            mp += dmp
            err = self.calcMass(mp) - m

        # bisection method between mpLast and mp
        mpa = mpLast
        mpb = mp
        mpc = (mpa + mpb) / 2

        tol = 1e-9
        errc = 2*tol

        while abs(errc) > tol:

            erra = self.calcMass(mpa) - m
            errb = self.calcMass(mpb) - m
            errc = self.calcMass(mpc) - m

            if (erra*errc) < 0:
                mpb = mpc

            elif (errb*errc) < 0:
                mpa = mpc

            else:
                print("Bisection method has failed, no crossing was present")
                exit()

            mpc = (mpa + mpb) / 2

        self.mp = mpc
        self.pFrac = self.mp / m
        self.mf = m - self.mp
        self.burnTime = self.mp / self.mDot


    def update(self, params:dict) -> None:

        if params.get('total_mass') is not None:
            self.mtot = params['total_mass']
            self.solveTotalMass()

        if params.get('isp') is not None:
            self.isp = params['isp']

        if params.get('mass_flow') is not None:
            self.mDot = params['mass_flow']

        self.thrust = self.isp * 9.80665 * self.mDot



class Rocket():
    # a rocket consists of one or more stages. These stages are defined by their propellant and inert masses, as well as their propulsion system's performance
    # the rocket is used to carry a payload of some mass, and provides functions for calculating vehicle performance combined with a payload

    def __init__(self, stages:list['Stage']) -> None:
        self.stages = stages
        self.canVary, self.freeStages = self.getVariableStages()
        self.currentStageNum = 0
        
        

    def getMass(self) -> float:

        m = 0

        for stage in self.stages:
            m += stage.mtot

        return m

    
    def getDv(self, payload:float) -> list['float']:

        dv = np.zeros((len(self.stages)), float)
        m0 = payload

        for i in range(len(self.stages) - 1, - 1, -1):
            
            m0 += self.stages[i].mtot
            mf = m0 - self.stages[i].mp
            dv[i] = self.stages[i].isp * 9.80665 * log(m0 / mf)

        return dv.tolist()
                

    def getAcceleration(self, payload:float, n:int) -> tuple['float', 'float']:
        """Returns the minimum and maximum accelerations for stage index n"""

        m0 = payload

        for i in range(n, len(self.stages)):
            m0 += self.stages[i].mtot
            
        aMin = self.stages[n].thrust / m0
        aMax = self.stages[n].thrust / (m0 - self.stages[n].mp)

        return aMin, aMax
    

    def getVariableStages(self):

        variableStages = []
        freeStages = []

        for i in range(0, len(self.stages)):
            variableStages.append(self.stages[i].canVary)
            
            if self.stages[i].canVary:
                freeStages.append(i)

        return variableStages, freeStages
    

    def updateStage(self, n:int, params:dict) -> None:
        self.stages[n].update(params)

    
    def scoreStages(self, payload:float) -> bool:

        self.stageScores = []
        total = 0.0

        for i in range(0, len(self.stages)):

            stage = self.stages[i]
            score = {}

            massLimits = stage.constraints.get('mass')
            accLimits = stage.constraints.get('acceleration')
            timeLimits = stage.constraints.get('time')
            
            if massLimits is not None:
                score['mass'] = max(max(0, massLimits[0] - stage.mtot), max(0, stage.mtot - massLimits[1]))
                total += score['mass']

            if accLimits is not None:
                minAcc, maxAcc = self.getAcceleration(payload, i)
                score['acceleration'] = max(max(0, accLimits[0] - minAcc), max(0, maxAcc - accLimits[1]))
                total += score['acceleration']
            
            if timeLimits is not None:
                score['time'] = max(max(0, timeLimits[0] - stage.burnTime), max(0, stage.burnTime - timeLimits[1]))
                total += score['time']

            self.stageScores.append(score)

        return total == 0.0
    

    def getCurrentStage(self) -> Stage:
        return self.stages[self.currentStageNum]


    def advanceStage(self):
        self.currentStageNum = min(len(self.stages) - 1, self.currentStageNum + 1) # cannot advance beyond the final stage

      
    def getThrust(self, throttle:float) -> float:
        stg = self.getCurrentStage()
        return stg.thrust*throttle if stg.mp > 0 else 0
     

    def simulate(self, dt:float, isBurning:bool, throttle:float) -> None:
        """Allows the vehicles state to be updated over one timestep"""
        stg = self.getCurrentStage()
        stg.simulate(dt, isBurning, throttle)