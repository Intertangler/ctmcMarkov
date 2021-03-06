import matplotlib.pyplot as plt
import random
import numpy as np


stateLabels = np.arange(14)
def transitionRates(rate_on,rate_off,rate_mono_to_double,rate_double_to_mono):
    transitionRateMatrix = np.zeros((len(stateLabels),len(stateLabels)))
    # go from empty state to monovalent, single occupancy
    transitionRateMatrix[0, 1] = rate_on # molecules /site/sec, depends on soln conc.
    transitionRateMatrix[0, 2] = rate_on
    transitionRateMatrix[0, 3] = 0
    # go from monovalent, single occupancy state to empty state
    transitionRateMatrix[1, 0] = rate_off # molecules /site/sec, binding mechnx dependnt
    transitionRateMatrix[2, 0] = rate_off
    transitionRateMatrix[3, 0] = 0
    # go from monovalent, single occupancy to monovalent, double occupancy state
    transitionRateMatrix[1, 4] = rate_on # similar to all monovalent trans
    transitionRateMatrix[1, 6] = 0
    transitionRateMatrix[2, 4] = rate_on
    transitionRateMatrix[2, 5] = 0
    transitionRateMatrix[3, 5] = 0
    transitionRateMatrix[3, 6] = 0
    # reverse of above block
    transitionRateMatrix[4, 1] = rate_off 
    transitionRateMatrix[6, 1] = 0
    transitionRateMatrix[4, 2] = rate_off
    transitionRateMatrix[5, 2] = 0
    transitionRateMatrix[5, 3] = 0
    transitionRateMatrix[6, 3] = 0
    
    # go from monovalent, single occupancy to bivalent 
    transitionRateMatrix[1, 8] = rate_mono_to_double # interested in solving for this internal trans
    transitionRateMatrix[1, 10] = 0
    transitionRateMatrix[2, 8] = rate_mono_to_double
    transitionRateMatrix[2, 9] = 0
    transitionRateMatrix[3, 9] = 0
    transitionRateMatrix[3, 10] = 0
    # reverse of above block
    transitionRateMatrix[8, 1] = rate_double_to_mono 
    transitionRateMatrix[10, 1] = 0
    transitionRateMatrix[8, 2] = rate_double_to_mono
    transitionRateMatrix[9, 2] = 0
    transitionRateMatrix[9, 3] = 0
    transitionRateMatrix[10, 3] = 0
    
    # go from monovalent, double occupancy to monovalent triple occupancy
    transitionRateMatrix[4, 7] = 0 # molecules /site/sec, depends on soln conc.
    transitionRateMatrix[5, 7] = 0
    transitionRateMatrix[6, 7] = 0
    # reverse of above block
    transitionRateMatrix[7, 4] = 0 # molecules /site/sec, binding mechnx dependnt
    transitionRateMatrix[7, 5] = 0
    transitionRateMatrix[7, 6] = 0
    
    # go from monovalent, double occupancy to bivalent/monovalent
    transitionRateMatrix[4, 12] = 0
    transitionRateMatrix[4, 13] = 0
    transitionRateMatrix[5, 11] = 0
    transitionRateMatrix[5, 13] = 0
    transitionRateMatrix[6, 12] = 0
    transitionRateMatrix[6, 11] = 0
    # reverse of above block
    transitionRateMatrix[12, 4] = 0
    transitionRateMatrix[13, 4] = 0
    transitionRateMatrix[11, 5] = 0
    transitionRateMatrix[13, 5] = 0
    transitionRateMatrix[12, 6] = 0
    transitionRateMatrix[11, 6] = 0
    
    # go from bivalent to bivalent/monovalent
    transitionRateMatrix[8, 11] = 0 # molecules /site/sec, depends on soln conc.
    transitionRateMatrix[9, 12] = 0
    transitionRateMatrix[10, 13] = 0
    # reverse of above block
    transitionRateMatrix[11, 8] = 0 # molecules /site/sec, binding mechnx dependnt
    transitionRateMatrix[12, 9] = 0
    transitionRateMatrix[13, 10] = 0
    return transitionRateMatrix

def exitRates(transitionRateMatrix):
    exitRateArray = np.zeros((len(transitionRateMatrix)))
    for i in range(0,len(exitRateArray)):
        exitRateArray[i] = np.sum(transitionRateMatrix[i,:])
    return exitRateArray
    
def embeddedDTMCfn(transitionRateMatrix):
    embeddedDTMC = np.zeros((len(transitionRateMatrix),len(transitionRateMatrix)))
    for i in range(0,len(embeddedDTMC[0,:])):
        for j in range(0,len(embeddedDTMC[:,0])):
            if exitRateArray[j] == 0:
                embeddedDTMC[i,j] = 0
            else:
                embeddedDTMC[i,j] = transitionRateMatrix[j,i]/exitRateArray[j]
    return embeddedDTMC

def cumEmbeddedDTMCfn(transitionRateMatrix,embeddedDTMC):
    cumEmbeddedDTMC = np.zeros((len(transitionRateMatrix),len(transitionRateMatrix)))
    for i in range(0,len(embeddedDTMC[0,:])):
        for j in range(0,len(embeddedDTMC[:,0])):
            cumEmbeddedDTMC[j,i] = np.sum(embeddedDTMC[0:j+1,i])
            if 1 in cumEmbeddedDTMC[0:j,i]:
                cumEmbeddedDTMC[j,i]=0
            else:
                pass
            if cumEmbeddedDTMC[j,i] in cumEmbeddedDTMC[0:j,i]:
                cumEmbeddedDTMC[j,i]=0
    return cumEmbeddedDTMC



numberOfSites = 1000
timeSamples = 500
timestep_size = .5
samplingInterval = timeSamples*timestep_size
timepoint1 = samplingInterval/4


for run in range(0,10):
    k_on = 0.0001*10**run/2
    k_off = .1
    k_mono_bi = .2
    k_bi_mono = .01
    occupancyRecord = np.zeros((timeSamples))
    for site in range(0, numberOfSites):
        cumulative_walking_time = 0
        step = 0
        print site
        currentState = 0
        transitionRateMatrix = transitionRates(k_on,k_off,k_mono_bi,k_bi_mono)
        exitRateArray = exitRates(transitionRateMatrix)
        embeddedDTMC = embeddedDTMCfn(transitionRateMatrix)
        cumEmbeddedDTMC = cumEmbeddedDTMCfn(transitionRateMatrix,embeddedDTMC)
        walkRecord = np.zeros((5000,3))
        timepoint = samplingInterval/4
        while cumulative_walking_time < timepoint:
            randomP_exit = random.uniform(0,1)
            if exitRateArray[currentState] == 0:
                walkRecord[step][0] = currentState
                holdingTime = timepoint-cumulative_walking_time
                walkRecord[step][1] = holdingTime
                cumulative_walking_time = timepoint
                walkRecord[step][2] = cumulative_walking_time
                step += 1
                break
            else:
                pass
            holdingTime = np.log(1/(1-randomP_exit))/exitRateArray[currentState]
            randomP_partition = random.uniform(0,1)
            walkRecord[step][0] = currentState
            walkRecord[step][1] = holdingTime
            cumulative_walking_time = cumulative_walking_time + holdingTime
            if cumulative_walking_time > timepoint:
                walkRecord[step][2] = timepoint
                walkRecord[step][1] = timepoint - walkRecord[step-1][2]
                step += 1
                break
            else:
                walkRecord[step][2] = cumulative_walking_time
            for i in range(0, len(cumEmbeddedDTMC[:,0])):
                if randomP_partition < cumEmbeddedDTMC[i,currentState] and cumEmbeddedDTMC[i,currentState] != 0:
                    currentState = i
                    step += 1
                    break
                else:
                    pass






        transitionRateMatrix = transitionRates(0, k_off, k_mono_bi, k_bi_mono)
        exitRateArray = exitRates(transitionRateMatrix)
        embeddedDTMC = embeddedDTMCfn(transitionRateMatrix)
        cumEmbeddedDTMC = cumEmbeddedDTMCfn(transitionRateMatrix,embeddedDTMC)
        timepoint = samplingInterval
        while cumulative_walking_time < timepoint:
            randomP_exit = random.uniform(0,1)
            if exitRateArray[currentState] == 0:
                walkRecord[step][0] = currentState
                holdingTime = timepoint-cumulative_walking_time
                walkRecord[step][1] = holdingTime
                cumulative_walking_time = timepoint
                walkRecord[step][2] = cumulative_walking_time
                step += 1
                break
            else:
                pass
            holdingTime = np.log(1/(1-randomP_exit))/exitRateArray[currentState]
            randomP_partition = random.uniform(0,1)
            walkRecord[step][0] = currentState
            walkRecord[step][1] = holdingTime
            cumulative_walking_time = cumulative_walking_time + holdingTime
            if cumulative_walking_time > timepoint:
                walkRecord[step][2] = timepoint
                walkRecord[step][1] = timepoint - walkRecord[step-1][2]
                step += 1
                break
            else:
                walkRecord[step][2] = cumulative_walking_time
            for i in range(0, len(cumEmbeddedDTMC[:,0])):
                if randomP_partition < cumEmbeddedDTMC[i,currentState] and cumEmbeddedDTMC[i,currentState] != 0:
                    currentState = i
                    step += 1
                    break
                else:
                    pass



        # rec all the occupancy events for this particular site now that the records have been built in the last two loops
        timeNow = 0
        currentOccupancy = 0
        for t_step in range(1, timeSamples):
            timeNow = t_step*timestep_size
            lookupstate = walkRecord[:,0][np.argmin((walkRecord[:,2]) < timeNow)]
            if lookupstate == 0:
                currentOccupancy = 0
            elif lookupstate == 1 or lookupstate == 2 or lookupstate == 3 or lookupstate == 8 or lookupstate == 9 or lookupstate == 10:
                currentOccupancy = 1
            elif lookupstate == 7:
                currentOccupancy = 3
            else:
                currentOccupancy = 2
            occupancyRecord[t_step] = occupancyRecord[t_step] + currentOccupancy
    #        print occupancyRecord[step]

    norm_occupancyRecord = occupancyRecord/numberOfSites
    timeRange = np.arange(0,samplingInterval,timestep_size)
    plt.plot(timeRange,norm_occupancyRecord,lw=2,label="$k_{on} \, =$ " + str(k_on) + '$v \, mol^{-1}  t^{-1}$' + "$k_{off} \, =$ " + str(k_off) + ' $  t^{-1}$, ' + "$k_{mono2bi} \, =$ " + str(k_mono_bi) + ' $  t^{-1}$, ' + "$k_{bi2mono} \, =$ " + str(k_bi_mono) + ' $  t^{-1}$, ') # vol~{-1} \, mol^{-1}
    plt.legend()
    plt.rcParams.update({'font.size': 22})
    plt.xlabel("time")
    plt.ylabel("Mean Occupancy per Site")
    plt.ylim((0,3))
plt.show("k_on_.01_to_100"+".png")
plt.close()

