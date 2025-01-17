import PhysicsConstants as PC
from SimulationDataset import *
import numpy as np
import random


# Times in seconds
def ActivityAtTime( StartingActivity, TimeElapsed, HalfLife ):
    return StartingActivity * ( 2.0 ** ( -TimeElapsed / HalfLife ) )

def F18ActivityAtTime( StartingActivity, TimeElapsed ):
    return ActivityAtTime( StartingActivity, TimeElapsed, 109.77*60.0 )

def Zr89ActivityAtTime( StartingActivity, TimeElapsed ):
    return ActivityAtTime( StartingActivity, TimeElapsed, 78.41*60.0*60.0 )

def TracerActivityAtTime( StartingActivity, TimeElapsed, Isotope ):
    return StartingActivity * ( 2.0 ** ( -TimeElapsed / PC.halflives[Isotope] ) )



# Simulate poisson-distributed random decay times
# Based on example https://timeseriesreasoning.com/contents/poisson-process/
def DeltaT( DecayRate ):
    randUniform = random.random()
    return -math.log( 1.0 - randUniform ) / DecayRate

def ActivityTimeline( DecayRate, EndTime ):
    decayTimes = []
    time = 0.0
    while time < EndTime:
        time += DeltaT( DecayRate )
        if time < EndTime:
            decayTimes.append( time )

    return decayTimes

# Avoid pre-generating timelines since it blows up the RAM
def GenerateCoincidences( DecayRates, EndTime, TimeWindow, TimelinesOut = None ):
    coincidences = []
    coincidenceTimes = []

    # Set the first event in each timeline
    nextTimes = []
    for channel in range( len( DecayRates ) ):
        nextTimes.append( DeltaT( DecayRates[ channel ] ) )
        if TimelinesOut is not None:
            TimelinesOut.append( [] )

    unfinishedTimeline = True
    while unfinishedTimeline:

        # Check if we have anything to process
        unfinishedTimeline = False
        for nextTime in nextTimes:
            if nextTime <= EndTime:
                unfinishedTimeline = True
                break

        # Start the window at the earliest event
        time = min( nextTimes )

        # Find all events in the window
        coincidence = []
        for channelIndex, nextTime in enumerate( nextTimes ):

            while nextTime <= time + TimeWindow and nextTime <= EndTime:

                # Store the time point
                coincidence.append( channelIndex )
                if TimelinesOut is not None:
                    TimelinesOut[ channelIndex ].append( nextTime )

                # Update to next time point
                nextTime += DeltaT( DecayRates[ channelIndex ] )

            # Store the first event outside the window, for the next coincidence
            nextTimes[ channelIndex ] = nextTime

        # The last entry is empty, because all times now past end
        if len( coincidence ) > 0:
            coincidences.append( coincidence )
            if TimelinesOut is not None:
                coincidenceTimes.append( time )
    return coincidences, coincidenceTimes




# Closer to the NEMA calculation, hopefully similar results
def NECRFromHistogram( bins, values, SimulationWindow ):
    # Assume one bin edge is on the 20mm boundary
    lowValue = 0.0
    lowBin = 0
    highValue = 0.0
    highBin = 0
    histogramTotal = 0.0
    centralTotal = 0.0
    for binIndex in range( len( bins ) - 1 ):
        if bins[binIndex] == -20.0:
            lowValue = ( values[binIndex] + values[binIndex-1] ) / 2.0
            lowBin = binIndex
        elif bins[binIndex] == 20.0:
            highValue = ( values[binIndex] + values[binIndex-1] ) / 2.0
            highBin = binIndex
    
        if bins[binIndex] >= -120.0 and bins[binIndex] < 120.0:
            histogramTotal += values[binIndex]
        if bins[binIndex] >= -20.0 and bins[binIndex] < 20.0:
            centralTotal += values[binIndex]
            
    centralBackground = ( lowValue + highValue ) / 2.0
    centralBackground *= ( highBin - lowBin ) # number of bins in the region
    
    trues = centralTotal - centralBackground
    rPlusS = histogramTotal - trues
    necr = 0.0
    if histogramTotal > 0.0:
        necr = trues * trues / histogramTotal

    return necr/SimulationWindow, trues/SimulationWindow, rPlusS/SimulationWindow

import matplotlib.pyplot as mpl

# base the coincidence window on whether the decay actually triggers the detector
def DetectedCoincidences( DecayRates, DecayData, SimulationWindow, CoincidenceWindow, DetectorRadius, ZMin=0.0, ZMax=0.0,
                          UsePhotonTime=False, EnergyResolution=0.0, TimeResolution=0.0 ):
    hitRadii = []
    sameEventIDArray = []
    trueEvents = 0
    allEvents = 0
    eventsOutsideMid = 0
    totalPhotons = 0
    usedPhotons = 0
    recycledPhotons = 0
    windowsFromRecycled = 0
    recyclePhotons = []

    # Set the first event in each timeline
    nextTimes = []
    for channel in range( len( DecayRates ) ):
        nextTimes.append( DeltaT( DecayRates[ channel ] ) )

    unfinishedTimeline = True
    while unfinishedTimeline:

        # Check if we have anything to process
        unfinishedTimeline = False
        for nextTime in nextTimes:
            if nextTime <= SimulationWindow:
                unfinishedTimeline = True
                break
        #Start the window at the earliest event
        time = min( nextTimes )
        event = [] # Have to start empty, then add. Otherwise end up taking a reference and bad things happen

        # Fill the event with first photon
        if UsePhotonTime and len(recyclePhotons) > 0:
            recyclePhotons = sorted( recyclePhotons )
            firstRecycleTime, firstRecyclePhoton = recyclePhotons[0]
            if firstRecycleTime < time:
                time = firstRecycleTime
                event += [firstRecyclePhoton]
                recyclePhotons.remove((firstRecycleTime, firstRecyclePhoton))
                windowsFromRecycled += 1

        # Without a recycled photon, generate a new one
        if len(event) == 0:
            minChannel = np.argmin( nextTimes )
            event += DecayData[ minChannel ].SampleOneEvent( EnergyResolution, TimeResolution )
            nextTimes[ minChannel ] += DeltaT( DecayRates[ minChannel ] )

        # Calculate the times the photons actually reach the detector
        detectionTimes = []
        if UsePhotonTime:
            for photon in event:
                detectionTimes.append( time + (photon[DATASET_TIME]*1E-9) )

        # Don't start a coincidence window if the event didn't pass cuts
        if len( event ) == 0:
            continue

        # Build event from all decays in the window
        for channelIndex, nextTime in enumerate( nextTimes ):

            while nextTime <= time + CoincidenceWindow and nextTime <= SimulationWindow:

                # Store the time point
                newDecay = DecayData[ channelIndex ].SampleOneEvent( EnergyResolution, TimeResolution )
                event += newDecay
                if UsePhotonTime:
                    for photon in newDecay:
                        detectionTimes.append( nextTime + (photon[DATASET_TIME]*1E-9) )

                # Update to next time point
                nextTime += DeltaT( DecayRates[ channelIndex ] )

            # Store the first event outside the window, for the next coincidence
            nextTimes[ channelIndex ] = nextTime

        # The last entry is empty, because all times now past end
        if len( event ) == 0:
            continue

        # Trim photons that reach the detector too late
        if UsePhotonTime:
            firstDetection = min( detectionTimes )
            trimmedEvent = []
            for i, photon in enumerate( event ):
                totalPhotons += 1
                if detectionTimes[i] < firstDetection + CoincidenceWindow:
                    trimmedEvent += [photon]
                    usedPhotons += 1
                else:
                    recyclePhotons.append((detectionTimes[i], photon))

            #Try recycling
            for detectionTime, photon in recyclePhotons:
                if detectionTime < time:
                    # We have passed this photon in the timeline, so get rid of it
                    # In theory this shouldn't happen
                    recyclePhotons.remove((detectionTime, photon))
                elif detectionTime >= time and detectionTime <= time + CoincidenceWindow:
                    trimmedEvent += [photon]
                    recyclePhotons.remove((detectionTime, photon))
                    recycledPhotons += 1

            event = trimmedEvent

        # Classify the events
        if TwoHitEvent( event, DetectorRadius, ZMin, ZMax ):
            hitRadius = FindHitRadius( event, DetectorRadius )
            hitRadii.append( hitRadius )
            hasSameEventID = SameEventID( event )
            sameEventIDArray.append( hasSameEventID )

            #check if event is a true, scatter or random
             
            #allEvents += 1
            #if BackToBackEvent( event, DetectorRadius, ZMin, ZMax ):
            #    trueEvents += 1
            #else:
            #    eventsOutsideMid += 1
    # Simplistic NECR calculation
    # The true events are those in the inner radius (2cm) minus expected background
    # background expectation is within the outer radius (12cm)
    # so true is (inside 2cm) - (inside 12cm, outside 2)/5 assuming a uniform background
    # not full NEMA calculation but comparable
    #trueEvents -= ( float(eventsOutsideMid) / 5.0 )
    #randomPlusScatter = allEvents - trueEvents
    #necr = 0
    #if allEvents > 0:
    #    necr = trueEvents * trueEvents / allEvents
    #
    # Return the hit radii for a better NECR calculation later
    #return necr/SimulationWindow, trueEvents/SimulationWindow, randomPlusScatter/SimulationWindow, hitRadii

    #if UsePhotonTime:
    #    lostPhotons = totalPhotons - usedPhotons
    #    print( str( lostPhotons ) + " photons lost of " + str( totalPhotons ) + " total (" + str( lostPhotons * 100 / totalPhotons ) + "%) and " + str(recycledPhotons) + " recycled, " + str(windowsFromRecycled) + " recycled to new window" )

    y, x, patches = mpl.hist( hitRadii, bins=26, range=[-130, 130] )
    mpl.clf()
    # First find the true and R+S coincidences
    NECRInSimWin, truesInSimWin, rPlusSInSimWin = NECRFromHistogram( x, y, SimulationWindow)
    # Then, find the number of scattered events by re-running NECRFromHistograms only for photons that have the same EventID
    # Find indexes of events in which both photons have the same eventID
    hitRadiiCoinc = []
    for i, val in enumerate(sameEventIDArray):
        if val == True:
            hitRadiiCoinc.append(hitRadii[i])
    y, x, patches = mpl.hist( hitRadiiCoinc, bins=26, range=[-130, 130] )
    mpl.clf()
    NECRTmp, truesTmp, scattersInSimWin = NECRFromHistogram( x, y, SimulationWindow)
    randomsInSimWin = rPlusSInSimWin - scattersInSimWin

    return NECRInSimWin, truesInSimWin, rPlusSInSimWin, scattersInSimWin, randomsInSimWin

# A generator to be used in different analysis methods
# will yield events in sequence
def GenerateCoincidences_new( DecayRates, DecayData, SimulationWindow, CoincidenceWindow, UsePhotonTime=False, EnergyResolution=0.0, TimeResolution=0.0 ):

    recyclePhotons = []

    # Set the first event in each timeline
    nextTimes = []
    for channel in range( len( DecayRates ) ):
        nextTimes.append( DeltaT( DecayRates[ channel ] ) )

    unfinishedTimeline = True
    while unfinishedTimeline:

        #Start the window at the earliest event
        windowStartTime = min( nextTimes )
        event = [] # Have to start empty, then add. Otherwise end up taking a reference and bad things happen
        detectionTimes = []

        # Fill the event with first photon
        if UsePhotonTime and len(recyclePhotons) > 0:
            recyclePhotons = sorted( recyclePhotons )
            firstRecycleTime, firstRecyclePhoton = recyclePhotons[0]
            if firstRecycleTime < windowStartTime:
                recyclePhotons.remove((firstRecycleTime, firstRecyclePhoton))
                windowStartTime = firstRecycleTime
                firstRecyclePhoton[DATASET_TIME] = firstRecycleTime
                detectionTimes.append( firstRecycleTime )
                event += [firstRecyclePhoton]

        # Without a recycled photon, generate new ones
        if len(event) == 0:
            minChannel = np.argmin( nextTimes )
            event += DecayData[ minChannel ].SampleOneEvent( EnergyResolution, TimeResolution )
            nextTimes[ minChannel ] += DeltaT( DecayRates[ minChannel ] )

            # Calculate the times the photons actually reach the detector
            if UsePhotonTime:
                for i, photon in enumerate( event ):
                    photonTime = windowStartTime + (photon[DATASET_TIME]*1E-9)
                    detectionTimes.append( photonTime )
                    event[i][DATASET_TIME] = photonTime

        # Don't start a coincidence window if the event didn't pass cuts
        if len( event ) == 0:
            continue

        # Build event from all decays in the window
        for channelIndex, nextTime in enumerate( nextTimes ):

            while nextTime <= windowStartTime + CoincidenceWindow and nextTime <= SimulationWindow:

                # Store the time point
                newDecay = DecayData[ channelIndex ].SampleOneEvent( EnergyResolution, TimeResolution )
                if UsePhotonTime:
                    for i, photon in enumerate( newDecay ):
                        photonTime = nextTime + (photon[DATASET_TIME]*1E-9)
                        detectionTimes.append( photonTime )
                        newDecay[i][DATASET_TIME] = photonTime

                event += newDecay

                # Update to next time point
                nextTime += DeltaT( DecayRates[ channelIndex ] )

            # Store the first event outside the window, for the next coincidence
            nextTimes[ channelIndex ] = nextTime

        # Check if we have anything to process in future
        unfinishedTimeline = False
        for nextTime in nextTimes:
            if nextTime <= SimulationWindow:
                unfinishedTimeline = True
                break

        if UsePhotonTime:

            # Trim photons that reach the detector too late
            firstDetection = min( detectionTimes )
            trimmedEvent = []
            for i, photon in enumerate( event ):
                if detectionTimes[i] < firstDetection + CoincidenceWindow:
                    trimmedEvent += [photon]
                else:
                    recyclePhotons.append((detectionTimes[i], photon))

            # Try recycling photons into the event
            for detectionTime, photon in recyclePhotons:
                if detectionTime < windowStartTime:
                    # We have passed this photon in the timeline, so get rid of it
                    # In theory this shouldn't happen
                    print( "WARNING: recycled photon was never used" )
                    recyclePhotons.remove((detectionTime, photon))
                elif detectionTime <= windowStartTime + CoincidenceWindow:
                    trimmedEvent += [photon]
                    recyclePhotons.remove((detectionTime, photon))

            event = trimmedEvent

        # Avoid (unlikely?) case where all photons trimmed
        if len( event ) == 0:
            continue

        yield event

# base the coincidence window on whether the decay actually triggers the detector
def DetectedCoincidences_new( DecayRates, DecayData, SimulationWindow, CoincidenceWindow, DetectorRadius, ZMin=0.0, ZMax=0.0,
                          UsePhotonTime=False, EnergyResolution=0.0, TimeResolution=0.0 ):

    hitRadii = []
    sameEventIDArray = []

    # Use the new generator method (in the python sense) for coincidences
    for event in GenerateCoincidences_new( DecayRates, DecayData, SimulationWindow, CoincidenceWindow, UsePhotonTime, EnergyResolution, TimeResolution ):

        # Classify the events
        if TwoHitEvent( event, DetectorRadius, ZMin, ZMax ):
            hitRadius = FindHitRadius( event, DetectorRadius )
            hitRadii.append( hitRadius )
            hasSameEventID = SameEventID( event )
            sameEventIDArray.append( hasSameEventID )

    y, x, patches = mpl.hist( hitRadii, bins=26, range=[-130, 130] )
    mpl.clf()
    # First find the true and R+S coincidences
    NECRInSimWin, truesInSimWin, rPlusSInSimWin = NECRFromHistogram( x, y, SimulationWindow)
    # Then, find the number of scattered events by re-running NECRFromHistograms only for photons that have the same EventID
    # Find indexes of events in which both photons have the same eventID
    hitRadiiCoinc = []
    for i, val in enumerate(sameEventIDArray):
        if val == True:
            hitRadiiCoinc.append(hitRadii[i])
    y, x, patches = mpl.hist( hitRadiiCoinc, bins=26, range=[-130, 130] )
    mpl.clf()
    NECRTmp, truesTmp, scattersInSimWin = NECRFromHistogram( x, y, SimulationWindow)
    randomsInSimWin = rPlusSInSimWin - scattersInSimWin

    return NECRInSimWin, truesInSimWin, rPlusSInSimWin, scattersInSimWin, randomsInSimWin
