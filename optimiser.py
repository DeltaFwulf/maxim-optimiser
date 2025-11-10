import numpy as np
from math import cos, pi, sqrt
from scipy.optimize import minimize

from rocket import Rocket
from copy import deepcopy



def scoreRocket(angs:list['float'], rocket:Rocket, mass:float, payload:float) -> float:

    fixedMass = 0
    for stage in rocket.stages:
        fixedMass += stage.m0 if not stage.canVary else 0.0

    freeMass = mass - fixedMass


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


    stageMasses = angs2masses(totalMass=freeMass, angs=angs)

    for i in range(0, len(rocket.freeStages)):
        params = {'total_mass':stageMasses[i]}
        rocket.updateStage(n=rocket.freeStages[i], params=params)

    rocket.scoreStages(payload=payload)

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

    dv = np.sum(rocket.getDv(payload))
    return -dv + massPenalty + accPenalty + timePenalty


def optimiseStages(rocket:Rocket, mRocket:float, mPayload:float) -> Rocket:
    """Calculates a rocket with a given total mass that delivers the highest possible delta-v according to design constraints."""

    # can this rocket be optimised? if not, return the same rocket immediately
    if len(rocket.freeStages) <= 1:
        return rocket

    # set input limits:
    angLims = [-pi/2, pi/2]
    bounds = []


    for i in range(0, len(rocket.freeStages)-1):#
        bounds.append(angLims)

    bounds = tuple(bounds)
    
    x0 = (0.5 * np.ones((len(rocket.freeStages) - 1), float)).tolist()
  
    minimize(fun=scoreRocket, x0=x0, args=(rocket, mRocket, mPayload), method='Nelder-Mead', bounds=bounds, tol=1e-9)

    return rocket



def setMass(rocketIn:Rocket, rocketMass:float, payload:float) -> Rocket:
    """Optimises a rocket's stages given a fixed vehicle and payload mass"""

    rocket = rocketIn
    mpl = payload
    
    if len(rocket.freeStages) > 1:
        optimiseStages(rocket, mRocket=rocketMass, mPayload=mpl)
        print(f"Rocket optimised")

    else:
        print("Rocket cannot be optimised, it only has one stage. Returning vehicle performance")
    

    return rocket



def lightestRocket(rocketIn:Rocket, payload:float, dv:float) -> Rocket:
    """Optimises the lightest possible rocket to carry a given payload with specified dv requirement"""

    mArr = np.logspace(base=10, start=2, stop=5, num=30, dtype=float)
    
    for i in range(mArr.size):

        # try different masses until something exceeds the dv requirement, then run shooting method between the two values:
        rocket = optimiseStages(deepcopy(rocketIn), mArr[i], payload)
        err = np.sum(rocket.getDv(payload)) - dv

        if err > 0:
            break

    m0 = mArr[i-1]
    m1 = mArr[i]

    tol = 1e-6

    err0 = 2*tol
    err1 = 2*tol

    while abs(err1) > tol:

        rocket0 = optimiseStages(deepcopy(rocketIn), m0, payload)
        rocket1 = optimiseStages(deepcopy(rocketIn), m1, payload)

        err0 = np.sum(rocket0.getDv(payload)) - dv
        err1 = np.sum(rocket1.getDv(payload)) - dv

        m2 = m0 + (err0 / (err0 - err1))*(m1 - m0)

        m0 = m1
        m1 = m2

    return rocket1


def maxPayload(rocketIn:Rocket, rocketMass:float, dv:float):
    """Determines the maximum payload that a rocket of given mass can carry with the specified dv. Returns the optimised rocket and its maximum payload."""

    p0 = 1
    p1 = 100

    tol = 1e-6

    err0 = 2*tol
    err1 = 2*tol

    rocketUnladen = optimiseStages(deepcopy(rocketIn), rocketMass, 0.0)
    if np.sum(rocketUnladen.getDv(0.0)) < dv:
        return rocketIn, 0.0

    while abs(err1) > tol:

        rocket0 = optimiseStages(deepcopy(rocketIn), rocketMass, p0)
        rocket1 = optimiseStages(deepcopy(rocketIn), rocketMass, p1)

        err0 = np.sum(rocket0.getDv(p0)) - dv
        err1 = np.sum(rocket1.getDv(p1)) - dv

        p2 = p0 + (err0 / (err0 - err1))*(p1 - p0)
        p0 = p1
        p1 = p2

    return rocket1, p1









# def maxPayload(rocket, mMin:float, mMax:float, dv:float):

#     # Calculate the maximum payload of a vehicle with given dv requirements against vehicle total scale
#     # for now, assume that the propellant mass fraction is held constant, but introduce a scaling law later

#     # NOTE: the rocket's mass excludes the payload mass, this is assumed to be on top of a standard vehicle design

#     # also, plot the curve using matplotlib to give a trend



#     m = np.linspace(mMin, mMax, 10)
#     maxPayload = np.zeros_like(m)

#     for i in range(m.size):

#         if len(rocket.stages) > 1: # optimise the stages for this mass
        
#             tol = 1e-9
#             err = [2*tol, 2*tol]
#             pl = [0, 100]

#             while abs(err) > tol:

#                 # optimise the payload such that an optimised multistage rocket can just deliver the required delta v when at the total mass.
#                 roc0 = deepcopy(rocket)
#                 roc1 = deepcopy(rocket)

#                 # FOR NOW: increase the mass of the rocket output by the payload returned... this is probably a bad idea though

#                 roc0.payload = pl[0]
#                 roc1.payload = pl[1]

#                 roc0 = targetDv(roc0, m[i])
#                 roc1 = targetDv(roc1, m[i])+

#                 err[0] = m[i] - roc0.getMass()
#                 err[1] = m[i] - roc1.getMass()

#                 pl2 = pl[0] + (err[0] / (err[0] - err[1]))*(pl[1] - pl[0])
#                 pl[0] = pl[1]
#                 pl[1] = pl2

#             maxPayload = pl2

#         else:
#             # direct calculation
#             pFrac = exp(-dv / ())