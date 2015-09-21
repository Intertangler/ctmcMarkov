__author__ = 'Group Administrator'
import networkx as nx
import matplotlib.pyplot as plt
import random
import csv
import numpy as np
import pygraphviz as pvg

#G=pvg.AGraph(strict=False,directed=True)
#nodelist=[
#            'empty_00',
#            'monovalent_single_occupancy_001',
#            'monovalent_single_occupancy_002',
#            'monovalent_single_occupancy_003',
#            'monovalent_double_occupancy_004',
#            'monovalent_double_occupancy_005',
#            'monovalent_double_occupancy_006',
#            'monovalent_triple_occupancy_007',
#            'bivalent_008',
#            'bivalent_09',
#            'bivalent_10',
#            'bivalent_monovalent_11',
#            'bivalent_monovalent_12',
#            'bivalent_monovalent_13',
#            ]
#G.add_nodes_from(nodelist)
#G.node_attr['shape']='box'
#G.node_attr['color']='#D0D0D0'
#G.node_attr['style']='filled'
#G.add_edge('empty_01','monovalent_single_occupancy_002',color='black')
#G.add_edge('monovalent_single_occupancy_002','empty_01',color='black')
#G.layout()
#G.layout(prog='fdp')
#G.draw('3_sites.png')









rate_on = 2 # molecules per site per second
rate_off = 1 # per second
rate_mono_to_double = 2 # binding events per second
rate_double_to_mono = 1 # unbinding events per second

stateLabels = np.arange(14)
def transitionRates(rate_on,rate_off,rate_mono_to_double,rate_double_to_mono):
    transitionRateMatrix = np.zeros((len(stateLabels),len(stateLabels)))
    # go from empty state to monovalent, single occupancy
    transitionRateMatrix[0, 1] = rate_on*3 # molecules /site/sec, depends on soln conc.
    transitionRateMatrix[0, 2] = rate_on*3
    transitionRateMatrix[0, 3] = rate_on*3
    # go from monovalent, single occupancy state to empty state
    transitionRateMatrix[1, 0] = rate_off # molecules /site/sec, binding mechnx dependnt
    transitionRateMatrix[2, 0] = rate_off
    transitionRateMatrix[3, 0] = rate_off
    # go from monovalent, single occupancy to monovalent, double occupancy state
    transitionRateMatrix[1, 4] = rate_on*2 # similar to all monovalent trans
    transitionRateMatrix[1, 6] = rate_on*2
    transitionRateMatrix[2, 4] = rate_on*2
    transitionRateMatrix[2, 5] = rate_on*2
    transitionRateMatrix[3, 5] = rate_on*2
    transitionRateMatrix[3, 6] = rate_on*2
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

def get_index(arr, val):                                                                
        index = np.searchsorted(arr, val)                                                            
        if arr[index] == val:                                                                        
            return index 
            
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
# initialize in unbound state
timeSteps = 1000
stepSize = .001
timeline = np.arange(0,timeSteps*stepSize,stepSize)
occupancyRecord = np.zeros((timeSteps))
totalSites = 1000

for site in range(0,totalSites):
    print site
    currentState = 0
    
    clocks = np.zeros((len(transitionRateMatrix)))           
    mainClock = 0
    iterations = 500
    walkRecord = np.zeros((2*iterations-1,3))
    
    # random walk
    transitionRateMatrix = transitionRates(2,1,2,1)
    exitRateArray = exitRates(transitionRateMatrix)
    embeddedDTMC = embeddedDTMCfn(transitionRateMatrix)
    cumEmbeddedDTMC = cumEmbeddedDTMCfn(transitionRateMatrix,embeddedDTMC)
    for step in range(1,iterations):
        randomP_exit = random.uniform(0,1)
        holdingTime = np.log(1/(1-randomP_exit))/exitRateArray[currentState]
        randomP_partition = random.uniform(0,1)
        for i in range(0, len(cumEmbeddedDTMC[:,0])):
            if randomP_partition < cumEmbeddedDTMC[i,currentState] and cumEmbeddedDTMC[i,currentState] != 0:
                walkRecord[step][0] = currentState
                walkRecord[step][1] = holdingTime
                clocks[currentState] = clocks[currentState] + holdingTime
                mainClock = mainClock + holdingTime
                walkRecord[step][2] = mainClock
#                print holdingTime
                currentState = i
                break
            else:
                pass
    # random walk 2 
    transitionRateMatrix = transitionRates(0,1,2,1)
    exitRateArray = exitRates(transitionRateMatrix)
    embeddedDTMC = embeddedDTMCfn(transitionRateMatrix)
    cumEmbeddedDTMC = cumEmbeddedDTMCfn(transitionRateMatrix,embeddedDTMC)
    for step in range(iterations+1,2*iterations):
        randomP_exit = random.uniform(0,1)
        holdingTime = np.log(1/(1-randomP_exit))/exitRateArray[currentState]
        randomP_partition = random.uniform(0,1)
        for i in range(0, len(cumEmbeddedDTMC[:,0])):
            if randomP_partition < cumEmbeddedDTMC[i,currentState] and cumEmbeddedDTMC[i,currentState] != 0:
                walkRecord[step][0] = currentState
                walkRecord[step][1] = holdingTime
                clocks[currentState] = clocks[currentState] + holdingTime
                mainClock = mainClock + holdingTime
                walkRecord[step][2] = mainClock
#                print holdingTime
                currentState = i
                break
            else:
                pass
    # discrete time averaging - can we get the exact occupancy? or at least truncated
    # analogous to spr which samples detected amount each discrete time interval
    
        
    timeNow = 0
    currentOccupancy = 0
    for step in range(1, timeSteps):
        timeNow = step*stepSize
        currentWalkerTime = find_nearest(walkRecord[:,2], timeNow)
        currentState = walkRecord[:,0][(get_index(walkRecord[:,2],currentWalkerTime))]
        if currentState == 0:
            currentOccupancy = 0
        elif currentState == 1 or currentState == 2 or currentState == 3 or currentState == 8 or currentState == 9 or currentState == 10:
            currentOccupancy = 1
        elif currentState == 7:
            currentOccupancy = 3
        else:
            currentOccupancy = 2
        occupancyRecord[step] = occupancyRecord[step] + currentOccupancy
#        print occupancyRecord[step]






plt.figure()
plt.plot(timeline,occupancyRecord)
plt.show()



#plt.figure()
#partitionSpace = plt.bar(range(len(transitionRateMatrix)),clocks,
#                 alpha=0.4,
#                 color='black')
#plt.show()
            
## divide up a range from 0 to 1 according to transition probabilities so that a random number will fall into a category
#partitionedSpaceMatrix = transitionProbabilityMatrix
#for i in range(len(transitionProbabilityMatrix)):
#    for j in range(len(transitionProbabilityMatrix[i, :])):
#        if j == 0:
#            partitionedSpaceMatrix[i, j] = 0 + transitionProbabilityMatrix[i, j]
#        else:
#            partitionedSpaceMatrix[i, j] = partitionedSpaceMatrix[i, j-1] + transitionProbabilityMatrix[i, j]
#
## choose a random number and execute transition according to current state and ...
## available neighbors - this is analogous to picking the particle's energy ...
## should this number be changed to sample from the boltzmann distribution with temperature as a parameter as well?
#
#counts = np.zeros(len(partitionedSpaceMatrix[0, :]))
#
## we want to know the minimum value which the randomvariable is less than in the partition array ...
## by subtracting the random variable from the upper limit of each bin, we get either a negative number or positive ...
## if we throw out the negatives and find the minimum positive number, that is the bin we want to sort into
#for i in range(0, 100000, 1):
#    randomVariable = random.random()
#    scoreArray = partitionedSpaceMatrix[currentState, :] - randomVariable
#    for j in range(len(scoreArray)):
#        if scoreArray[j] <= 0:
#            scoreArray[j] = 1
#        else:
#            pass
#    currentState = np.argmin(scoreArray)
#    counts[currentState] = counts[currentState] + 1
#    normalizedCounts = counts/np.max(counts)
#
#
