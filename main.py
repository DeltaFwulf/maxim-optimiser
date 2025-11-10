from rocket import Stage, Rocket
from optimiser import setMass, lightestRocket, maxPayload
import numpy as np
import matplotlib.pyplot as plt

# TODO: create propellant class, allow mixes to be defined to give propellant bulk density 
# TODO: if rocket has only one free stage, the total mass should be corrected to the enforced mass

def main():

    mRocket = 20000
    dvMin = 5000
    dvMax = 20000

    stdFixedMass = 100.0
    stgLambda = [0.9, 0.9, 0.9, 0.9, 0.9]
    stgIsp = [325, 325, 325, 325, 325]

    # define the rockets you want to compare here
    stagesA = []
    stagesA.append(Stage(mTotal=mRocket, propFrac=stgLambda[0], isp=stgIsp[0], mDot=10, canVary=True, fixedMass=stdFixedMass, constraints={}))

    stagesB = []
    stagesB.append(Stage(mTotal=0.0, propFrac=stgLambda[0], isp=stgIsp[0], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesB.append(Stage(mTotal=0.0, propFrac=stgLambda[1], isp=stgIsp[1], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))

    stagesC = []
    stagesC.append(Stage(mTotal=0.0, propFrac=stgLambda[0], isp=stgIsp[0], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesC.append(Stage(mTotal=0.0, propFrac=stgLambda[1], isp=stgIsp[1], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesC.append(Stage(mTotal=0.0, propFrac=stgLambda[2], isp=stgIsp[2], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))

    stagesD = []
    stagesD.append(Stage(mTotal=0.0, propFrac=stgLambda[0], isp=stgIsp[0], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesD.append(Stage(mTotal=0.0, propFrac=stgLambda[1], isp=stgIsp[1], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesD.append(Stage(mTotal=0.0, propFrac=stgLambda[2], isp=stgIsp[2], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesD.append(Stage(mTotal=0.0, propFrac=stgLambda[3], isp=stgIsp[3], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))

    stagesE = []
    stagesE.append(Stage(mTotal=0.0, propFrac=stgLambda[0], isp=stgIsp[0], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesE.append(Stage(mTotal=0.0, propFrac=stgLambda[1], isp=stgIsp[1], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesE.append(Stage(mTotal=0.0, propFrac=stgLambda[2], isp=stgIsp[2], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesE.append(Stage(mTotal=0.0, propFrac=stgLambda[3], isp=stgIsp[3], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))
    stagesE.append(Stage(mTotal=0.0, propFrac=stgLambda[4], isp=stgIsp[4], mDot=55.0, canVary=True, fixedMass=stdFixedMass, constraints={}))

    rocketA = Rocket(stages=stagesA)
    rocketB = Rocket(stages=stagesB)
    rocketC = Rocket(stages=stagesC)
    rocketD = Rocket(stages=stagesD)
    rocketE = Rocket(stages=stagesE)

    rocketA.name = 'one stage'
    rocketB.name = 'two stages'
    rocketC.name = 'three stages'
    rocketD.name = 'four stages'
    rocketE.name = 'five stages'

    rockets = {rocketA.name:rocketA, rocketB.name:rocketB, rocketC.name:rocketC, rocketD.name:rocketD, rocketE.name:rocketE}
    maxPayloadVsDv(rockets, mRocket, [dvMin, dvMax])



def optimiseGivenMass():

    stages = []
    # stages.append(Stage(mTotal=4342, propFrac=0.92, isp=296, mDot=50, canVary=False, fixedMass=50, constraints={'mass':[1000, 1e8]}))
    # stages.append(Stage(mTotal=950, propFrac=0.87, isp=287, mDot=11, canVary=False, fixedMass=50, constraints={'mass':[200, 1e8]}))
    # stages.append(Stage(mTotal=936, propFrac=0.7, isp=380, mDot=5, canVary=True, fixedMass=50, constraints={'mass':[200, 1e8]}))

    stages.append(Stage(mTotal=2000.0, propFrac=0.8, isp=300.0, mDot=10.0, canVary=True, fixedMass=0.0, constraints={}))
    stages.append(Stage(mTotal=200.0, propFrac=0.8, isp=300.0, mDot=5.0, canVary=True, fixedMass=0.0, constraints={}))

    rocket = Rocket(stages=stages)
    rocketMass = 2000
    payload = 100.0

    rocket = setMass(rocket, rocketMass, payload)

    for j in range(0, len(rocket.stages)):
        print(f"Stage {j}:\n=========================================================")
        print(f"Total Mass:\t\t{'%.3f' % rocket.stages[j].m0} kg")
        print(f"Propellant Mass:\t{'%.3f' % rocket.stages[j].mp} kg")
        print(f"Final Mass:\t\t{'%.3f' % rocket.stages[j].mf} kg")
        print(f"Required burn time:\t{'%.3f' % rocket.stages[j].burnTime} s")
        print(f"delta-v:\t\t{'%.3f' % rocket.getDv(payload=payload)[j]} m/s")
        print('\n')

    print(f"Total dv: {'%.3f' % np.sum(rocket.getDv(payload))}")



def minViableMass():

    stages = []
    stages.append(Stage(mTotal=4000, propFrac=0.92, isp=296, mDot=50, canVary=True, fixedMass=50, constraints={'mass':[4292, 4292], 'acceleration':[20, 60]}))
    stages.append(Stage(mTotal=900, propFrac=0.876, isp=287, mDot=11, canVary=True, fixedMass=50, constraints={'acceleration':[20, 60]}))
    stages.append(Stage(mTotal=250, propFrac=0.7, isp=380, mDot=5, canVary=True, fixedMass=50, constraints={'acceleration':[20, 60]}))
   
    rocket = Rocket(stages=stages)
    payload = 10
    dv = 8000

    rocket = lightestRocket(rocket, payload, dv)
    
    print("Rocket optimised")
    print(f"Minimum mass: {rocket.getMass()} kg")

    for j in range(0, len(rocket.stages)):
        print(f"Stage {j}:\n=========================================================")
        print(f"Total Mass:          {'%.3f' % rocket.stages[j].m0} kg")
        print(f"Propellant Mass:     {'%.3f' % rocket.stages[j].mp} kg")
        print(f"Final Mass:          {'%.3f' % rocket.stages[j].mf} kg")
        print(f"Required burn time:  {'%.3f' % rocket.stages[j].burnTime} s")
        print(f"delta-v:             {'%.3f' % rocket.getDv(payload=payload)[j]} m/s")
        print(f"Minimum accleration: {'%.3f' % (rocket.getAcceleration(payload, j)[0])} m/s^2")
        print(f"Maximum accleration: {'%.3f' % (rocket.getAcceleration(payload, j)[1])} m/s^2")
        print('\n')



def maxPayloadVsMass():

    stages = []

    stages.append(Stage(mTotal=2000.0, propFrac=0.8, isp=300.0, mDot=10.0, canVary=True, fixedMass=0.0, constraints={}))
    stages.append(Stage(mTotal=200.0, propFrac=0.8, isp=300.0, mDot=10.0, canVary=True, fixedMass=0.0, constraints={}))

    rocket = Rocket(stages=stages)

    dv = 5000

    mArr = np.linspace(2000, 3000, 100)
    mpl = np.zeros_like(mArr)
    
    for i in range(mArr.size):
        rocketOut, mpl[i] = maxPayload(rocket, mArr[i], dv)


    fig, ax = plt.subplots()
    ax.plot(mArr, mpl, '-k')
    ax.set_xlabel("Rocket mass, kg")
    ax.set_ylabel("Max Payload, kg")
    ax.set_title(f"Maximum payload vs rocket mass for dv = {'%.1f' % dv} m/s")
    ax.grid(True)

    plt.show()



def maxPayloadVsDv(rockets:list['Rocket'], mRocket:float, dvRange:list['float']) -> None:

    dvArr = np.linspace(dvRange[0], dvRange[1], 50)
    payloads = {}

    for key in rockets:

        rocket = rockets[key]

        pl = np.zeros_like(dvArr)
        for i in range(dvArr.size):
            _, pl[i] = maxPayload(rocket, mRocket, dvArr[i])

        payloads.update({key:pl})

    # how many stages is best?
    bestNum = np.zeros_like(dvArr)
    for i in range(dvArr.size):
        
        maxPl = 0.0
        for key in payloads:
            
            attempt = payloads[key][i]

            if attempt > maxPl:
                maxPl = attempt     
                best = len(rockets.get(key).stages)

        bestNum[i] = best

    fig, ax = plt.subplots()
    for key in payloads:
        payloads[key] = payloads[key][payloads[key] > 0]
        ax.plot(dvArr[:payloads[key].size], payloads[key], '-', label=key)

    ax.set_xlabel("dv, m/s")
    ax.set_ylabel("Max Payload, kg")
    ax.set_title(f"Maximum payload vs dv for rocket with total mass {'%.1f' % mRocket} kg")
    ax.legend()
    ax.grid(True)

    fig2, ax2 = plt.subplots()
    ax2.plot(dvArr, bestNum, '-k', label="best # stages")
    ax2.set_xlabel("dv, m/s")
    ax2.set_ylabel("Optimum stage number")
    ax2.set_title("Optimum number of stages against rocket dv requirement")
    
    plt.show()



main()