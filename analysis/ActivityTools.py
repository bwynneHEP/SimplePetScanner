# Times in seconds
def ActivityAtTime( StartingActivity, TimeElapsed, HalfLife ):
    return StartingActivity * ( 2.0 ** ( -TimeElapsed / HalfLife ) )

def F18ActivityAtTime( StartingActivity, TimeElapsed ):
    return ActivityAtTime( StartingActivity, TimeElapsed, 109.77*60.0 )

def Zr89ActivityAtTime( StartingActivity, TimeElapsed ):
    return ActivityAtTime( StartingActivity, TimeElapsed, 78.41*60.0*60.0 )


import random

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


from SimulationDataset import *
import numpy as np

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
    
    true = centralTotal - centralBackground
    rPlusS = histogramTotal - true
    necr = true * true / histogramTotal
    
    return necr/SimulationWindow, true/SimulationWindow, rPlusS/SimulationWindow

import matplotlib.pyplot as mpl

# base the coincidence window on whether the decay actually triggers the detector
def DetectedCoincidences( DecayRates, DecayData, SimulationWindow, CoincidenceWindow, DetectorRadius, ZMin=0.0, ZMax=0.0 ):

    hitRadii = []
    trueEvents = 0
    allEvents = 0
    eventsOutsideMid = 0

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

        # Start the window at the earliest event
        time = min( nextTimes )
        minChannel = np.argmin( nextTimes )
        event = [] # Have to start empty, then add. Otherwise end up taking a reference and bad things happen
        event += DecayData[ minChannel ].SampleOneEvent()
        nextTimes[ minChannel ] += DeltaT( DecayRates[ minChannel ] )

        # TODO investigate the event time itself (column index 2)
        # probably doesn't matter, we're just trying to use random decays

        # Don't start a coincidence window if the event didn't pass cuts
        if len( event ) == 0:
          continue

        # Build event from all decays in the window
        for channelIndex, nextTime in enumerate( nextTimes ):

            while nextTime <= time + CoincidenceWindow and nextTime <= SimulationWindow:

                # Store the time point
                event += DecayData[ channelIndex ].SampleOneEvent()

                # Update to next time point
                nextTime += DeltaT( DecayRates[ channelIndex ] )

            # Store the first event outside the window, for the next coincidence
            nextTimes[ channelIndex ] = nextTime

        # The last entry is empty, because all times now past end
        if len( event ) == 0:
            continue

        # Classify the events
        if TwoHitEvent( event, DetectorRadius, ZMin, ZMax ):
            hitRadii.append( FindHitRadius( event, DetectorRadius ) )
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

    y, x, patches = mpl.hist( hitRadii, bins=26, range=[-130, 130] )
    return NECRFromHistogram( x, y, SimulationWindow )

