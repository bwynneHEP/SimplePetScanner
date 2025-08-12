# A class for using the data in more-or-less the same format as it is loaded

import numpy as np
from SimulationDataset import DATASET_EVENT, DATASET_ENERGY, DATASET_TIME, DATASET_R, DATASET_PHI, DATASET_Z, DATASET_PHOTON_LENGTH

class LegacyDatasetReader:

  def __init__( self, InputDataset ):
    self.inputDataset = InputDataset
    self.unusedEvents = []
    self.usedEvents = []
    self.energyMin = InputDataset.energyMin
    self.energyMax = InputDataset.energyMax
    self.totalDecays = InputDataset.totalDecays

    # Allow for decays that weren't detected
    for i in range( self.totalDecays ):
      self.unusedEvents.append( i )


  def ReferenceOneEvent( self, RNG=None ):

    # NOTE: ideally this is an internal method that just returns event data by reference

    # Check if we have any events left
    if len( self.unusedEvents ) == 0:
      self.unusedEvents = self.usedEvents
      if RNG == None:
        self.inputDataset.RNG.shuffle( self.unusedEvents )
      else:
        RNG.shuffle( self.unusedEvents )
      self.usedEvents = []

    eventID = self.unusedEvents.pop(-1)
    self.usedEvents.append( eventID )
    if eventID in self.inputDataset.inputData:
      return self.inputDataset.inputData[ eventID ]
    else:
      return []


  def SampleOneEvent( self, EnergyResolution=0.0, TimeResolution=0.0 ):

    # NOTE: older-style method for retrieving events one-at-a-time with resolution effects applied

    modifiedEvent = []
    for photon in self.ReferenceOneEvent():

      # Apply resolution effects to each measured photon
      newPhoton = [ value for value in photon ]
      newPhoton[DATASET_ENERGY] = photon[DATASET_ENERGY] * RNG.normal( 1.0, EnergyResolution ) # Energy resolution as a percentage
      newPhoton[DATASET_TIME] = photon[DATASET_TIME] + RNG.normal( 0.0, TimeResolution ) # Time resolution as absolute ns

      # Apply energy cut to modified photon
      keepPhoton = True
      if self.energyMin is not None and newPhoton[DATASET_ENERGY] < self.energyMin:
        keepPhoton = False
      if self.energyMax is not None and newPhoton[DATASET_ENERGY] > self.energyMax:
        keepPhoton = False
      if keepPhoton:
        modifiedEvent.append( newPhoton )

    return modifiedEvent


  def SampleEventsAtTimes( self, Times, RNG=None ):

    # NOTE: method for (hopefully faster) batch processing

    result = []

    for time in Times:
      for photon in self.ReferenceOneEvent( RNG ):

        newPhoton = [ value for value in photon ]
        newPhoton[DATASET_TIME] += ( time * 1e9 ) # convert to ns

        # Flattened across events
        result.append( newPhoton )

    return np.array( result )


  def size( self ):
    return self.totalDecays

  def GetModuleIDs( self ):
    return self.inputDataset.moduleMap
