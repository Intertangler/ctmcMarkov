import matplotlib.pyplot as plt
import random
import numpy as np


stateLabels = np.arange(14)
def transitionRates(rate_on,rate_off,rate_mono_to_double,rate_double_to_mono):
    transitionRateMatrix = np.zeros((len(stateLabels),len(stateLabels)))
    # go from empty state to monovalent, single occupancy
    transitionRateMatrix[0, 1] = rate_on # molecules /site/sec, depends on soln conc.
    transitionRateMatrix[0, 2] = rate_on
    transitionRateMatrix[0, 3] = rate_on
    # go from monovalent, single occupancy state to empty state
    transitionRateMatrix[1, 0] = rate_off # molecules /site/sec, binding mechnx dependnt
    transitionRateMatrix[2, 0] = rate_off
    transitionRateMatrix[3, 0] = rate_off
    # go from monovalent, single occupancy to monovalent, double occupancy state
    transitionRateMatrix[1, 4] = rate_on # similar to all monovalent trans
    transitionRateMatrix[1, 6] = rate_on
    transitionRateMatrix[2, 4] = rate_on
    transitionRateMatrix[2, 5] = rate_on
    transitionRateMatrix[3, 5] = rate_on
    transitionRateMatrix[3, 6] = rate_on
    # reverse of above block
    transitionRateMatrix[4, 1] = rate_off 
    transitionRateMatrix[6, 1] = rate_off
    transitionRateMatrix[4, 2] = rate_off
    transitionRateMatrix[5, 2] = rate_off
    transitionRateMatrix[5, 3] = rate_off
    transitionRateMatrix[6, 3] = rate_off
    
    # go from monovalent, single occupancy to bivalent 
    transitionRateMatrix[1, 8] = rate_mono_to_double # interested in solving for this internal trans
    transitionRateMatrix[1, 10] = rate_mono_to_double
    transitionRateMatrix[2, 8] = rate_mono_to_double
    transitionRateMatrix[2, 9] = rate_mono_to_double
    transitionRateMatrix[3, 9] = rate_mono_to_double
    transitionRateMatrix[3, 10] = rate_mono_to_double
    # reverse of above block
    transitionRateMatrix[8, 1] = rate_double_to_mono 
    transitionRateMatrix[10, 1] = rate_double_to_mono
    transitionRateMatrix[8, 2] = rate_double_to_mono
    transitionRateMatrix[9, 2] = rate_double_to_mono
    transitionRateMatrix[9, 3] = rate_double_to_mono
    transitionRateMatrix[10, 3] = rate_double_to_mono
    
    # go from monovalent, double occupancy to monovalent triple occupancy
    transitionRateMatrix[4, 7] = rate_on # molecules /site/sec, depends on soln conc.
    transitionRateMatrix[5, 7] = rate_on
    transitionRateMatrix[6, 7] = rate_on
    # reverse of above block
    transitionRateMatrix[7, 4] = rate_off # molecules /site/sec, binding mechnx dependnt
    transitionRateMatrix[7, 5] = rate_off
    transitionRateMatrix[7, 6] = rate_off
    
    # go from monovalent, double occupancy to bivalent/monovalent
    transitionRateMatrix[4, 12] = rate_mono_to_double 
    transitionRateMatrix[4, 13] = rate_mono_to_double
    transitionRateMatrix[5, 11] = rate_mono_to_double
    transitionRateMatrix[5, 13] = rate_mono_to_double
    transitionRateMatrix[6, 12] = rate_mono_to_double
    transitionRateMatrix[6, 11] = rate_mono_to_double
    # reverse of above block
    transitionRateMatrix[12, 4] = rate_double_to_mono 
    transitionRateMatrix[13, 4] = rate_double_to_mono
    transitionRateMatrix[11, 5] = rate_double_to_mono
    transitionRateMatrix[13, 5] = rate_double_to_mono
    transitionRateMatrix[12, 6] = rate_double_to_mono
    transitionRateMatrix[11, 6] = rate_double_to_mono
    
    # go from bivalent to bivalent/monovalent
    transitionRateMatrix[8, 11] = rate_on # molecules /site/sec, depends on soln conc.
    transitionRateMatrix[9, 12] = rate_on
    transitionRateMatrix[10, 13] = rate_on
    # reverse of above block
    transitionRateMatrix[11, 8] = rate_off # molecules /site/sec, binding mechnx dependnt
    transitionRateMatrix[12, 9] = rate_off
    transitionRateMatrix[13, 10] = rate_off
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




numberOfSites = 300
timeSamples = 500
timestep = 10
samplingInterval = timeSamples*timestep



# construct timeline of events up to the latest sampling time point and note the final state at that point
# for site in range(0,numberOfSites-1):

for run in range(0,1):
    k_on = 1
    occupancyRecord = np.zeros((timeSamples))
    for site in range(0, numberOfSites):
        cumulative_walking_time = 0
        step = 0
        print site
        currentState = 0
        transitionRateMatrix = transitionRates(.1,.1,.1,.1)
        exitRateArray = exitRates(transitionRateMatrix)
        embeddedDTMC = embeddedDTMCfn(transitionRateMatrix)
        cumEmbeddedDTMC = cumEmbeddedDTMCfn(transitionRateMatrix,embeddedDTMC)
        walkRecord = np.zeros((5000,3))
        while cumulative_walking_time < samplingInterval/4:
            randomP_exit = random.uniform(0,1)
            if exitRateArray[currentState] == 0:
#                holdingTime = samplingInterval/4-cumulative_walking_time
##                walkRecord[step][0] = currentState
#                walkRecord[step][1] = holdingTime
#                cumulative_walking_time += holdingTime
#                walkRecord[step-1][2] = cumulative_walking_time
#                currentState = currentState
#                step += 1
                break
            else:
                pass
            holdingTime = np.log(1/(1-randomP_exit))/exitRateArray[currentState]
            walkRecord[step][1] = holdingTime
            walkRecord[step][0] = currentState
            randomP_partition = random.uniform(0,1)
            for nextstate in range(0, len(cumEmbeddedDTMC[:,0])):
                if randomP_partition < cumEmbeddedDTMC[nextstate,currentState] and cumEmbeddedDTMC[nextstate,currentState] != 0:
                    # clocks[currentState] = clocks[currentState] + holdingTime
                    cumulative_walking_time += holdingTime
                    walkRecord[step][2] = cumulative_walking_time
                    currentState = nextstate
                    step += 1
                    break
                else:
                    pass
            
        cumulative_walking_time = samplingInterval/4
        

        print 'asdf' + str(currentState)
        transitionRateMatrix = transitionRates(0, .1, .1, .1)
        exitRateArray = exitRates(transitionRateMatrix)
        embeddedDTMC = embeddedDTMCfn(transitionRateMatrix)
        cumEmbeddedDTMC = cumEmbeddedDTMCfn(transitionRateMatrix,embeddedDTMC)
        while cumulative_walking_time < samplingInterval:
            randomP_exit = random.uniform(0,1)
            if exitRateArray[currentState] == 0:
                holdingTime = samplingInterval
                walkRecord[step][0] = currentState
                walkRecord[step][1] = holdingTime
                cumulative_walking_time += holdingTime
                walkRecord[step][2] = cumulative_walking_time
                step += 1
                break
            else:
                pass
            holdingTime = np.log(1/(1-randomP_exit))/exitRateArray[currentState]
            walkRecord[step][1] = holdingTime
            walkRecord[step][0] = currentState
            randomP_partition = random.uniform(0,1)
            for nextstate in range(0, len(cumEmbeddedDTMC[:,0])):
                if randomP_partition < cumEmbeddedDTMC[nextstate,currentState] and cumEmbeddedDTMC[nextstate,currentState] != 0:
                    # clocks[currentState] = clocks[currentState] + holdingTime
                    cumulative_walking_time += holdingTime
                    walkRecord[step][2] = cumulative_walking_time
                    currentState = nextstate
                    step += 1
                    break
                else:
                    pass

        # rec all the occupancy events for this particular site now that the records have been built in the last two loops
        timeNow = 0
        currentOccupancy = 0
        for tick in range(1, timeSamples):
            timeNow = tick*timestep
#            currentWalkerTime = np.argmin((walkRecord[:,2])>timeNow)
            lookupstate = walkRecord[:,0][np.argmin((walkRecord[:,2])<timeNow)]
            if lookupstate == 4 or lookupstate == 5 or lookupstate == 6 or lookupstate == 11 or lookupstate == 12 or lookupstate == 13:
                currentOccupancy = 2
            elif lookupstate == 1 or lookupstate == 2 or lookupstate == 3 or lookupstate == 8 or lookupstate == 9 or lookupstate == 10:
                currentOccupancy = 1
            elif lookupstate == 7:
                currentOccupancy = 3
            else:
                lookupstate = 0
            occupancyRecord[tick] = occupancyRecord[tick] + currentOccupancy/float(numberOfSites)
    #        print occupancyRecord[step]

#    norm_occupancyRecord = occupancyRecord/numberOfSites
    timeRange = np.arange(0,samplingInterval,timestep)
    plt.plot(timeRange,occupancyRecord)
plt.savefig("k_on_.01_to_100"+".png")
plt.close()

