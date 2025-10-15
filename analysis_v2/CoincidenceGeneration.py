from SimulationDataset import *

import numpy as np
#import time


# A method to create a set of time values corresponding to a given decay rate
# The values will be sufficient to cover the TimePeriod argument
# TimePeriod and DecayRate must have consistent units
# RNG is passed as an argument, following numpy recommendations
# There is effectively a "carry" operation in the form of StartTime:
#  this method will return the first event after the end of TimePeriod,
#  and this can then be passed as the first event in the next series
def TimeSeriesSingleChannel( TimePeriod, DecayRate, RNG, StartTime=0.0 ):

  # Multiply by reciprocal is still faster
  decayReciprocal = 1 / DecayRate

  # Guess how many events to generate
  eventsGuess = int( TimePeriod * DecayRate )
  if eventsGuess < 1:
    eventsGuess = 1

  # Entropy input
  randomArray = RNG.random( size=eventsGuess )

  # Convert to Poisson intervals
  unscaledDeltaT = -np.log( 1.0 - randomArray )

  # Scale by the expected decay rate
  scaledDeltaT = unscaledDeltaT * decayReciprocal
  
  # For batching, use the last event generated in the previous series as the first of this
  # The time will be presented as an "overflow" past the end of the last time period
  if StartTime > 0.0:
    scaledDeltaT = np.append( StartTime, scaledDeltaT )

  # Convert from delta-T to time series
  timeSeries = np.cumsum( scaledDeltaT )

  # For profiling
  #addonNumbers = []

  # Here's the unpleasant bit - make sure there are enough events to cover the time period
  while timeSeries[-1] < TimePeriod:
    
    # Guess how many to generate by progress so far
    progress = timeSeries[-1] / TimePeriod
    eventsGuess = int( len( timeSeries ) * ( 1.0 - progress ) / progress )

    # Set a minimum of 8 events - something something vectorisation
    # These cycles are more expensive than spare events
    if eventsGuess < 8:
      eventsGuess = 8
    #addonNumbers.append( eventsGuess )

    randomArrayExtra = RNG.random( size=eventsGuess )
    unscaledDeltaTExtra = -np.log( 1.0 - randomArrayExtra )
    scaledDeltaTExtra = unscaledDeltaTExtra * decayReciprocal
    
    # Offset from the end of existing series
    scaledDeltaTExtra[0] += timeSeries[-1]
    
    timeSeriesExtra = np.cumsum( scaledDeltaTExtra )

    # Append to existing series
    timeSeries = np.append( timeSeries, timeSeriesExtra )
 
  # Profiling results
  #print( "Addons: " + str( addonNumbers ) )
  #print( "Truncation: " + str( sum( timeSeries >= TimePeriod ) ) )

  # Work out where the first event in the next series should be
  startNext = timeSeries[ timeSeries >= TimePeriod ][0] - TimePeriod
  
  # Truncate to the time window (should always lose at least one event)
  timeSeries = timeSeries[ timeSeries < TimePeriod ]

  return timeSeries, startNext


# To generate multiple channels of decay data at different rates
# Specify a target (max) number of decays per channel with BatchSize
# The code will then calculate the time period that all decays
#  must populate, regardless of their decay rate
# This allows consistent combination of decay channels later
def TimeSeriesMultiChannel( BatchSize, DecayRates, RNG, StartTimes=None ):

  # Check inputs
  if not isinstance( DecayRates, np.ndarray ):
    DecayRates = np.array( DecayRates )
  if not DecayRates.ndim == 1:
    raise ValueError( "Decay rates must be a 1D float array with one entry per decay channel" )
  if not isinstance( StartTimes, np.ndarray ):
    StartTimes = np.zeros( len( DecayRates ) )
  assert len( StartTimes ) == len( DecayRates ), "One start time offset must be specified per channel"

  #start = time.time_ns()
  
  # It's not possible to have equal numbers of events in different channels
  # Even if they all had the same decay rate, the random element means they don't cover the same time range
  # So, define an expected time range from the fastest decay, and ensure all channels cover it
  fastestDecay = np.max( DecayRates )
  timeGuess = BatchSize / fastestDecay

  # Ragged array, one channel at a time
  # This is OK since the next operation is the per-channel photon gather anyway
  result = [[]] * len( DecayRates )
  for i in range( len( DecayRates ) ):
    result[i], startNext = TimeSeriesSingleChannel( timeGuess, DecayRates[i], RNG, StartTimes[i] )
    StartTimes[i] = startNext

  #end = time.time_ns()
  #print( "Make time series: " + str( end-start ) + "ns" )

  return result, timeGuess


# Take multiple decay channels and sample G4 photons for each channel
# Then merge all photons into a single timeline for coincidence calulation
def MergedPhotonStream( TimeSeries, DecayData, RNG, EnergyResolution=0.0, EnergyMin=0.0, EnergyMax=0.0, TimeResolution=0.0, Sort=True ):

  # Check inputs
  if len( TimeSeries ) != len( DecayData ):
    raise ValueError( "TimeSeries and DecayData must be arrays of the same (1st dimension) length" )

  #start = time.time_ns()

  # Convert decays into photons using data samples
  # Internally SampleEventsAtTimes() will flatten across decay events
  photons = DecayData[0].SampleEventsAtTimes( TimeSeries[0], RNG )

  # Flatten all decay channels into a single set of photons by appending
  for i in range( 1, len( TimeSeries ) ):

    morePhotons = DecayData[i].SampleEventsAtTimes( TimeSeries[i], RNG )

    if len( photons ) == 0:
      photons = morePhotons
    elif len( morePhotons ) > 0:
      photons = np.append( photons, morePhotons, axis=0 )

  # If there are no photons then skip the rest
  if len( photons ) == 0:
    return photons

  # Apply time resolution to photons
  if TimeResolution > 0.0:
    photonTimeOffsets = RNG.normal( 0, TimeResolution, size=len(photons) )
    photons[:,DATASET_TIME] += photonTimeOffsets

  # Apply energy resolution to photons
  if EnergyResolution > 0.0:
    photonEnergyMultipliers = RNG.normal( 1, EnergyResolution, size=len(photons) )
    photons[:,DATASET_ENERGY] *= photonEnergyMultipliers

  # Apply energy cut
  if EnergyMin > 0.0:
    photons = photons[ photons[:,DATASET_ENERGY] > EnergyMin ]
  if EnergyMax > 0.0:
    photons = photons[ photons[:,DATASET_ENERGY] < EnergyMax ]

  #Sort the photons by time
  if len( photons ) > 0 and Sort:
    photons = photons[ photons[:,DATASET_TIME].argsort() ]

  #end = time.time_ns()
  #print( "Create photon stream: " + str( end-start ) + "ns" )

  return photons


# Code that uses the above methods to generate a stream of coincidence windows
# In each window, details of all compatible photons are returned
# BatchSize allows control of the number of decays in memory at once,
#  but the code should run multiple batches until the SimulationWindow is
#  filled, i.e. the total time will be broken into an unknown number of
#  processing steps, each containing a (roughly) constant number of events
# DecayRates and DecayData describe each channel under consideration,
#  giving the rate of decays in that channel, and the dataset containing
#  the Geant4 data for those decays
# CoincidenceWindow is the time (ns) to collect photons for a coincidence,
#  and MultiWindow is a boolean flag permitting each photon to open a new
#  coincidence window (or not, if there is an existing window)
# TODO decide on time units - ns versus s (currently mixed)
def GenerateCoincidences( BatchSize, DecayRates, DecayData, RNG, CoincidenceWindow, SimulationWindow, MultiWindow, \
                          EnergyResolution=0.0, EnergyMin=0.0, EnergyMax=0.0, TimeResolution=0.0, ContinuousTimes=True ):

  # Since each decay is calculated by delta-T, use the last in each channel as an offset
  timeOffsets = np.zeros( len( DecayRates ) )

  # profiling
  #yieldCounter = 0
  #batchCounter = 0

  # Loop until end of the simulation
  totalTime = 0.0
  leftoverPhotons = None
  while totalTime < SimulationWindow:

    # Use the previous methods to create the photon timeline
    timeSeries, batchTimePeriod = TimeSeriesMultiChannel( BatchSize, DecayRates, RNG, timeOffsets )
    photonStream = MergedPhotonStream( timeSeries, DecayData, RNG, EnergyResolution, EnergyMin, EnergyMax, TimeResolution, Sort=False )

    #batchCounter += 1

    # Potentially this could just be an empty batch
    if len( photonStream ) == 0:
      continue

    # Photon times are generated from the start of the batch: convert to experiment time
    if ContinuousTimes:
      photonStream[:,DATASET_TIME] += totalTime

    # Re-use previous batch photons that were leftover
    if leftoverPhotons is not None and len( leftoverPhotons ) > 0:
      photonStream = np.append( leftoverPhotons, photonStream, axis=0 )

    # Sort photons after merging batches to avoid overlaps
    photonStream = photonStream[ photonStream[:,DATASET_TIME].argsort() ]

    # Find the last photon time, to make sure we don't overrun
    finalPhotonTime = photonStream[ -1, DATASET_TIME ]

    # Loop over the photons, opening coincidence windows
    startWindowIndex = 0
    while startWindowIndex < len( photonStream ):

      # Start the window with this photon
      thisPhoton = photonStream[startWindowIndex]
      thisPhotonTime = thisPhoton[DATASET_TIME]
      endWindowTime = thisPhotonTime + CoincidenceWindow

      # Truncate when the simulation is finished
      if ContinuousTimes and endWindowTime > SimulationWindow:
        #print( "Total batches: " + str( batchCounter ) )
        #print( "Yielded " + str( yieldCounter/batchCounter ) + " windows/batch" )
        return
      elif not ContinuousTimes and endWindowTime + totalTime > SimulationWindow:
        #print( "Total batches: " + str( batchCounter ) )
        #print( "Yielded " + str( yieldCounter/batchCounter ) + " windows/batch" )
        return

      # Check for needing a new batch
      # Original check was
      #if endWindowTime > finalPhotonTime:
      # However this fails if the new batch of photons includes entries that
      # should have fallen into an older window (due to TOF effects)
      # Adding at least a full window's distance prevents this
      if endWindowTime + CoincidenceWindow > finalPhotonTime:
        totalTime += ( batchTimePeriod * 1e9 )
        break

      # Create the window data by examining subsequent photons
      endWindowIndex = -1
      for nextPhotonIndex in range( startWindowIndex + 1, len( photonStream ) ):

        nextPhotonTime = photonStream[ nextPhotonIndex, DATASET_TIME ]

        # Check if photon is within the window
        if nextPhotonTime >= endWindowTime:
          endWindowIndex = nextPhotonIndex
          break
     
      #yieldCounter += 1

      yield photonStream[ startWindowIndex : endWindowIndex ]

      # Choose whether to allow multiple concurrent windows
      if MultiWindow:
        startWindowIndex += 1
      else:
        startWindowIndex = endWindowIndex

    # We've broken out of the loop because it's the end of a batch
    # Recycle the leftover photons for the next batch
    leftoverPhotons = photonStream[ startWindowIndex : ]

    # Offset times for a new batch
    if not ContinuousTimes:
      leftoverPhotons[ :, DATASET_TIME ] -= ( batchTimePeriod * 1e9 )
