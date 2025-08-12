# A class for (faster) data manipulation with zero-padded numpy arrays

import numpy as np
from SimulationDataset import DATASET_EVENT, DATASET_ENERGY, DATASET_TIME, DATASET_R, DATASET_PHI, DATASET_Z, DATASET_PHOTON_LENGTH

class NumpyDatasetReader:

  def __init__( self, InputDataset ):
    self.inputDataset = InputDataset
    self.unusedEvents = None
    self.energyMin = InputDataset.energyMin
    self.energyMax = InputDataset.energyMax
    self.totalDecays = InputDataset.totalDecays
    self.eventHitsMax = InputDataset.eventHitsMax
    self.sampleStartIndex = 0

    # Allow for decays that weren't detected
    self.unusedEvents = np.arange( self.totalDecays )

    # Make a numpy array to hold the data - alternative approach to sampling
    # NOTE: this is going to have a *lot* of zero-padding
    # Might not be worth the trade for the numpy speedup
    # Alternatively might need to use awkward arrays instead
    self.numpyData = np.zeros( [self.totalDecays, self.eventHitsMax, DATASET_PHOTON_LENGTH] )
    for eventIndex in self.inputDataset.inputData:
      for photonIndex in range( len( self.inputDataset.inputData[ eventIndex ] ) ):
        self.numpyData[ eventIndex, photonIndex, : ] = self.inputDataset.inputData[ eventIndex ][ photonIndex ]

    # Finished loading, so clear the input data
    self.inputDataset.inputData = None


  def SampleEventsAtTimes( self, Times, RNG=None ):

    # This first part of the code simply asks if the batch can be filled
    # If the batch is larger than the remaining events in the file,
    #  recursively calls this method twice and combines results
    # Note - won't work if the batch is larger than the whole file
    firstBatchMax = len( self.unusedEvents ) - self.sampleStartIndex
    if firstBatchMax < len( Times ):

      # Get as many events as you can from the current sequence
      firstPart = self.SampleEventsAtTimes( Times[:firstBatchMax] )

      # Reshuffle
      # Note it is a lot faster to shuffle this set of indices
      #  than to shuffle the actual data
      if RNG == None:
        self.RNG.shuffle( self.unusedEvents )
      else:
        RNG.shuffle( self.unusedEvents )
      self.sampleStartIndex = 0

      # Get remaining events
      secondPart = self.SampleEventsAtTimes( Times[firstBatchMax:] )

      # Combine the two sets of events
      return np.append( firstPart, secondPart, axis=0 )

    # Choose the event indices that will be used for the batch
    sampleIndices = self.unusedEvents[ self.sampleStartIndex : self.sampleStartIndex + len( Times ) ]
    self.sampleStartIndex = self.sampleStartIndex + len( Times )

    # Clone the event data corresponding to those indices
    events = self.numpyData[ sampleIndices ].copy()

    # Add the corresponding time offsets for each event
    events[ :, :, DATASET_TIME ] += ( Times * 1e9 )[ :, np.newaxis ] # convert to ns, then broadcast

    # Flatten across events to just give photons
    events = events.reshape( len( Times ) * self.eventHitsMax, DATASET_PHOTON_LENGTH )

    # Remove zero-padding
    events = events[ events[:,DATASET_ENERGY] > 0.0 ]
    return( events )


  def size( self ):
    return self.totalDecays

  def GetModuleIDs( self ):
    return self.inputDataset.moduleMap
