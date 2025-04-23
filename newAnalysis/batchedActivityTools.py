from SimulationDataset import *

import numpy as np
import time


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
  
  # Add the start time as an offset (allows batching)
  scaledDeltaT[0] += StartTime

  # Convert from delta-T to time series
  timeSeries = np.cumsum( scaledDeltaT )

  #addonNumbers = []

  # Here's the unpleasant bit - make sure there are enough events to cover the time period
  while timeSeries[-1] < ( TimePeriod + StartTime ):
    
    # Guess how many to generate by progress so far
    progress = ( timeSeries[-1] - StartTime ) / TimePeriod
    eventsGuess = int( len( timeSeries ) * ( 1.0 - progress ) )

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
 
  #print( "Addons: " + str( addonNumbers ) )
  #print( "Truncation: " + str( sum( timeSeries >= ( TimePeriod + StartTime ) ) ) )

  # Truncate to the time window (should always lose at least one event)
  # TODO check this
  #truncateFilter = timeSeries < ( TimePeriod + StartTime )
  #assert len( timeSeries ) - sum( truncateFilter ) > 1, "No events truncated: time series not filled"
  timeSeries = timeSeries[ timeSeries < ( TimePeriod + StartTime ) ]

  return timeSeries


def TimeSeriesMultiChannel( BatchSize, DecayRates, RNG, StartTime=0.0 ):

  # Check inputs
  if not isinstance( DecayRates, np.ndarray ):
    DecayRates = np.array( DecayRates )
  if not DecayRates.ndim == 1:
    raise ValueError( "Decay rates must be a 1D float array with one entry per decay channel" )
  
  start = time.time_ns()
  
  # It's not possible to have equal numbers of events in different channels
  # Even if they all had the same decay rate, the random element means they don't cover the same time range
  # So, define an expected time range from the fastest decay, and ensure all other channels cover it
  fastestDecay = np.max( DecayRates )
  timeGuess = BatchSize / fastestDecay

  # Ragged array, one channel at a time
  # This is OK since the next operation is the per-channel photon gather anyway
  result = [[]] * len( DecayRates )
  for i in range( len( DecayRates ) ):
    result[i] = TimeSeriesSingleChannel( timeGuess, DecayRates[i], RNG )

  end = time.time_ns()
  print( "Make time series: " + str( end-start ) + "ns" )

  return result


def MergedPhotonStream( TimeSeries, DecayData, RNG, EnergyResolution=0.0, EnergyMin=0.0, EnergyMax=0.0, TimeResolution=0.0 ):

  # Check inputs
  if len( TimeSeries ) != len( DecayData ):
    raise ValueError( "TimeSeries and DecayData must be arrays of the same (1st dimension) length" )
  for thing in DecayData:
    if not isinstance( thing, SimulationDataset ):
      raise ValueError( "DecayData must be instances of the SimulationDataset class" )

  start = time.time_ns()

  # Convert decays into photons using data samples
  photons = []
  for i, times in enumerate( TimeSeries ):
    # Uses += to flatten across channels
    # Internally will flatten across decay events
    photons += DecayData[i].SampleEventsAtTimes( times ) # TODO pass rng for photon sampling
  photons = np.array( photons )

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
  if len( photons ) > 0:
    photons = photons[ photons[:,DATASET_TIME].argsort() ]

  end = time.time_ns()
  print( "Create photon stream: " + str( end-start ) + "ns" )

  return photons


def GenerateCoincidences( BatchSize, DecayRates, DecayData, RNG, CoincidenceWindow, SimulationWindow, MultiWindow, \
                          EnergyResolution=0.0, EnergyMin=0.0, EnergyMax=0.0, TimeResolution=0.0 ):

  timeSeries = TimeSeriesMultiChannel( BatchSize, DecayRates, RNG, StartTime=0.0 )
  photonStream = MergedPhotonStream( timeSeries, DecayData, RNG, EnergyResolution, EnergyMin, EnergyMax, TimeResolution )

  # Find the last photon time, to make sure we don't overrun
  lastPhotonTime = photonStream[ -1, DATASET_TIME ]

  # Truncate when the simulation is finished
  if lastPhotonTime > SimulationWindow:
    lastPhotonTime = SimulationWindow
  endPhotonIndex = -1

  # Loop over the photons, opening coincidence windows
  startPhotonIndex = 0
  while startPhotonIndex < len( photonStream ):
    
    # Start the window with this photon
    thisPhoton = photonStream[startPhotonIndex]
    thisPhotonTime = thisPhoton[DATASET_TIME]
    endWindowTime = thisPhotonTime + CoincidenceWindow
    
    # Check for needing a new batch TODO: actually make new batches
    if endWindowTime > lastPhotonTime:
      endPhotonIndex = startPhotonIndex
      break

    # Create the window data by examining subsequent photons
    endWindowIndex = -1
    for nextPhotonIndex in range( startPhotonIndex + 1, len( photonStream ) ):

      nextPhotonTime = photonStream[ nextPhotonIndex, DATASET_TIME ]

      # Check if photon is within the window
      if nextPhotonTime >= endWindowTime:
        endWindowIndex = nextPhotonIndex
        break
    
    yield photonStream[ startPhotonIndex : endPhotonIndex ]

    # Choose whether to allow multiple concurrent windows
    if MultiWindow:
      startPhotonIndex += 1
    else:
      startPhotonIndex = endWindowIndex
