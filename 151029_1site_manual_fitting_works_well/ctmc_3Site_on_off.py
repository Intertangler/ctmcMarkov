import matplotlib.pyplot as plt
import random
import numpy as np
import time as time
from StringIO import StringIO
from scipy.optimize import minimize

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
    
def embeddedDTMCfn(transitionRateMatrix,exitRateArray):
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


def random_walk(concentration,step,timepoint,walkRecord, cumulative_walking_time,currentState):
        partition = np.zeros((len(stateLabels)))
        # print step
        transitionRateMatrix = transitionRates(k_on*concentration, k_off, k_mono_bi, k_bi_mono)
        exitRateArray = exitRates(transitionRateMatrix)
        embeddedDTMC = embeddedDTMCfn(transitionRateMatrix, exitRateArray)
        cumEmbeddedDTMC = cumEmbeddedDTMCfn(transitionRateMatrix,embeddedDTMC)
        while cumulative_walking_time < timepoint:
            randomP_exit = random.uniform(0,1)
            if exitRateArray[currentState] == 0:
                walkRecord[step][0] = currentState
                holdingTime = timepoint-cumulative_walking_time
                walkRecord[step][1] = holdingTime
                partition[currentState] = partition[currentState] + holdingTime
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
            partition[currentState] = partition[currentState] + holdingTime
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
        return walkRecord, currentState, step, cumulative_walking_time

def perform_series_run(rate_constants):
    with open('data.txt', "r") as dataFile:
        data = dataFile.read()
    actual_run = np.genfromtxt(StringIO(data), delimiter="\n")
    actual_run = actual_run/10
    numberOfSites = 50
    timeSamples = len(actual_run)
    timestep_size = 1
    samplingInterval = timeSamples*timestep_size
    for run in range(0,1):
        k_on = rate_constants[0]
        k_off = rate_constants[1]
        k_mono_bi = rate_constants[2]
        k_bi_mono = rate_constants[3]
        occupancyRecord = np.zeros((timeSamples))

        for site in range(0, numberOfSites):
            cumulative_walking_time = 0
            step = 0
            # print site
            currentState = 0
            walkRecord = np.zeros((5000, 3))
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=0, step=0, timepoint=84, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=0.025, step=step, timepoint=384, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=0, step=step, timepoint=475, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=0.05, step=step, timepoint=775, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=0, step=step, timepoint=866, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=.1, step=step, timepoint=1166, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=0, step=step, timepoint=1257, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=.25, step=step, timepoint=1557, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=0, step=step, timepoint=1656, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=.5, step=step, timepoint=1956, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)
            walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=0, step=step, timepoint=3000, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)

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
        sse = (norm_occupancyRecord-actual_run)**2
        master_error = np.sum(sse)
        if master_error < 500:
            plt.figure()
            plt.plot(timeRange,norm_occupancyRecord,lw=2,label="$k_{on} \, =$ " + str(k_on) + ' $v \, mol^{-1}  t^{-1}$, ' + "$k_{off} \, =$ " + str(k_off) + ' $  t^{-1}$, ' + "$k_{mono2bi} \, =$ " + str(k_mono_bi) + ' $  t^{-1}$, ' + "$k_{bi2mono} \, =$ " + str(k_bi_mono) + ' $  t^{-1}$, ') # vol~{-1} \, mol^{-1}

            plt.rcParams.update({'font.size': 7})
            plt.xlabel("time")
            plt.ylabel("Mean Occupancy per Site")
            plt.ylim((0,3))
            plt.plot(range(len(actual_run)),actual_run,lw=2,label='actual_run_data') # vol~{-1} \, mol^{-1}
            plt.legend()
            plt.savefig(str(time.time())+"occupancy_curve"+".png")
            plt.close()
        # plt.figure()
        # partitionSpace = plt.bar(range(len(stateLabels)),partition[0:14],
        #             alpha=0.4,
        #             color='black')
        # plt.savefig(str(time.time())+"partitionSpace"+".png")
        # plt.close()

    # plt.show("k_on_.01_to_100"+".png")
    # plt.close()



    return master_error



#intitial guess
k_on = .1
k_off = .05
k_mono_bi = 5
k_bi_mono = .05
rate_constants = np.array([k_on, k_off, k_mono_bi, k_bi_mono])
res = minimize(perform_series_run,rate_constants, method='Powell',
             options={'xtol': 1e-8, 'disp': True})
