10:23 2015-10-28
*beginning to look at how to fit models to alan's spr data
*first thing to adjust is the injection times - since his experiments are not simply on experiment and off experiment, but multiple concentrations
*best way may be to rewrite the main loop as a function so that i can call it with start times, end times, and concentrations, and walkrecord as inputs
**and output would be the modified walkrecord
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
					
12:41 2015-10-28
*function version of the random walk run works now
*inputs are concentration, step, timepoint, walkrecord, cumulative walking time, and the current state
*outputs are walk record, the current state, step, and the cumulative walking time
*to run the function, i set the return variables equal to the function with appropriate input parameters
    def random_walk(concentration,step,timepoint,walkRecord, cumulative_walking_time,currentState):
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


*running the function looks like:
	walkRecord, currentState, step, cumulative_walking_time = random_walk(concentration=15, step=0, timepoint=samplingInterval/2, walkRecord=walkRecord, cumulative_walking_time=cumulative_walking_time, currentState=currentState)

13:16 2015-10-28
*series of concentration steps is looking good now - ready to start working on the fitting part
