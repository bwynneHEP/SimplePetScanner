# File format defined in EnergyCounter.cpp:
#   std::cout << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " " << entry.first << " " << entry.second / keV << " ";
#   std::cout << m_averageTimeMap[ entry.first ] / ( entry.second * ns ) << " ";
#   std::cout << m_averageRMap[ entry.first ] / ( entry.second * mm ) << " ";
#   std::cout << m_averagePhiMap[ entry.first ] / entry.second << " ";
#   std::cout << m_averageZMap[ entry.first ] / ( entry.second * mm ) << std::endl;
# so
# EventID ModuleID Energy Time R Phi Z

# Avoid the module ID because it's never used and can be stringified
DATASET_EVENT, DATASET_ENERGY, DATASET_TIME, DATASET_R, DATASET_PHI, DATASET_Z, DATASET_PHOTON_LENGTH = 0, 1, 2, 3, 4, 5, 6


import math
import subprocess
import os
import numpy as np


class CsvFileReader:

  def __init__( self, InputPath ):
    self.inputFile = open( InputPath )
    self.nextLine = self.inputFile.readline()
    self.currentEvent = []

  def __del__( self ):
    self.inputFile.close()

  def Open( self ):
    return self.nextLine != ""

  def ReadHit( self ):
    if not self.Open():
      return None

    splitLine = self.nextLine.split(" ")

    # Assemble hit info
    eventID = int( splitLine[DATASET_EVENT] )
    wholeHit = [ eventID ] # Not floats
    for i in range( 2, len( splitLine ) ):
      wholeHit.append( float( splitLine[i] ) )

    self.nextLine = self.inputFile.readline()
    return wholeHit

  def ReadEvent( self ):
    if not self.Open():
      return None

    if len( self.currentEvent ) == 0:
      self.currentEvent.append( self.ReadHit() )

    while True:
      newHit = self.ReadHit()
      if newHit is None:
        return self.currentEvent

      if newHit[DATASET_EVENT] == self.currentEvent[0][DATASET_EVENT]:
        self.currentEvent.append( newHit )
      else:
        newEvent = self.currentEvent
        self.currentEvent = [newHit]
        return newEvent


class SimulationDataset:

  def __init__( self, InputPath, TotalDecays, EnergyMin=None, EnergyMax=None, ClusterLimitMM=None, RNG=None ):
    self.inputData = {}
    self.unusedEvents = None
    self.energyMin = EnergyMin
    self.energyMax = EnergyMax
    self.totalDecays = TotalDecays
    self.hitCount = 0
    self.eventHitsMax = 0
    self.sampleStartIndex = 0
    self.RNG = RNG
    if self.RNG == None:
      self.RNG = np.random.default_rng()

    if TotalDecays < 1:
      print( "ERROR: Requesting an empty dataset" )
      return

    currentEvent = -1
    eventCount = 0

    # Parse input
    inputFile = CsvFileReader( InputPath )
    while inputFile.Open():

      wholeHit = inputFile.ReadHit()
      eventID = wholeHit[ DATASET_EVENT ]

      # Multiple lines (hits) can go into a single event
      if eventID in self.inputData:
        self.AddHit( eventID, wholeHit, ClusterLimitMM )
      else:
        # Input file should be ordered and thus old event complete
        if currentEvent in self.inputData:
          self.eventHitsMax = max( self.eventHitsMax, len( self.inputData[currentEvent] ) )

        # Start a new event
        self.inputData[ eventID ] = [ wholeHit ]
        currentEvent = eventID
        eventCount += 1
        self.hitCount += 1

    print( str(eventCount) + " events loaded (" + str( self.totalDecays ) + " simulated) with average " + str( self.hitCount / self.totalDecays ) + " hits/event" )

    # Allow for decays that weren't detected
    self.unusedEvents = np.arange( self.totalDecays )

    # Make a numpy array to hold the data - alternative approach to sampling
    # NOTE: this is going to have a *lot* of zero-padding
    # Might not be worth the trade for the numpy speedup
    # Alternatively might need to use awkward arrays instead
    self.numpyData = np.zeros( [self.totalDecays, self.eventHitsMax, DATASET_PHOTON_LENGTH] )
    for eventIndex in self.inputData:
      for photonIndex in range( len( self.inputData[ eventIndex ] ) ):
        self.numpyData[ eventIndex, photonIndex, : ] = self.inputData[ eventIndex ][ photonIndex ]

    # Finished loading, so clear the input data
    self.inputData = None


  def AddHit( self, ExistingEventID, NewHit, ClusterLimitMM ):

    if ClusterLimitMM is None:
      self.inputData[ ExistingEventID ].append( NewHit )
      self.hitCount += 1

    else:
      oldEvent = self.inputData[ ExistingEventID ]
      newE = NewHit[DATASET_ENERGY]
      newZ = NewHit[DATASET_Z]
      newPhi = NewHit[DATASET_PHI]
      newR = NewHit[DATASET_R]
      keepHit = True

      # Attempt to add hit to each hit to the new event
      for oldHitIndex, oldHit in enumerate( oldEvent ):
        oldE = oldHit[DATASET_ENERGY]
        oldZ = oldHit[DATASET_Z]
        oldPhi = oldHit[DATASET_PHI]
        oldR = oldHit[DATASET_R]
        deltaZ = newZ-oldZ
        deltaPhiR = (newPhi*newR)-(oldPhi*oldR)
        delta = math.sqrt( deltaZ*deltaZ + deltaPhiR*deltaPhiR )

        # Merge hits within threshold
        if delta < ClusterLimitMM:
          mergedE = newE + oldE
          mergedZ = ( newE*newZ + oldE*oldZ ) / mergedE
          mergedPhi = ( newE*newPhi + oldE*oldPhi ) / mergedE
          mergedR = ( newE*newR + oldE*oldR ) / mergedE
          mergedT = ( newE*NewHit[DATASET_TIME] + oldE*oldHit[DATASET_TIME] ) / mergedE
          self.inputData[ ExistingEventID ][ oldHitIndex ] = [ oldHit[DATASET_EVENT], mergedE, mergedT, mergedR, mergedPhi, mergedZ ]
          keepHit = False
          break

      if keepHit:
        self.inputData[ ExistingEventID ].append( NewHit )
        self.hitCount += 1


  def SampleEventsAtTimes( self, Times, RNG=None ):

    #TODO check if time precision is sufficient. Make sure times don't get too big - only have to be relative to start of batch

    #TODO check this code for reshuffle
    firstBatchMax = len( self.unusedEvents ) - self.sampleStartIndex
    if firstBatchMax < len( Times ):

      # Get as many events as you can from the current sequence
      firstPart = self.SampleEventsAtTimes( Times[:firstBatchMax] )

      # Reshuffle
      if RNG == None:
        self.RNG.shuffle( self.unusedEvents )
      else:
        RNG.shuffle( self.unusedEvents )
      self.sampleStartIndex = 0

      # Get remaining events
      secondPart = self.SampleEventsAtTimes( Times[firstBatchMax:] )

      # Combine the two sets of events
      return np.append( firstPart, secondPart )

    sampleIndices = self.unusedEvents[ self.sampleStartIndex : self.sampleStartIndex + len( Times ) ]
    self.sampleStartIndex = self.sampleStartIndex + len( Times )
#    print( sampleIndices )
    events = self.numpyData[ sampleIndices ].copy()
#    print( events )
    events[ :, :, DATASET_TIME ] += ( Times * 1e9 )[ :, np.newaxis ] # convert to ns, then broadcast
#    print( events )
    events = events.reshape( len( Times ) * self.eventHitsMax, DATASET_PHOTON_LENGTH )
#    print( events )
    events = events[ events[:,DATASET_ENERGY] > 0.0 ]
#    print( events )
    return( events )


  def size( self ):
    return self.totalDecays


#
# End of the class, now just defining general methods
#

def SameEventID( Event ):
  return Event[0][DATASET_EVENT] == Event[1][DATASET_EVENT]


def FindHitRadius( Event, DetectorRadius ):

  if len( Event ) != 2:
    return -1.0

  # Calculate delta phi
  phi1 = Event[0][DATASET_PHI]
  phi2 = Event[1][DATASET_PHI]
  deltaPhi = phi1 - phi2
  while deltaPhi > math.pi:
    deltaPhi -= 2.0 * math.pi
  while deltaPhi < -math.pi:
    deltaPhi += 2.0 * math.pi

  if deltaPhi < 0.0:
    return -DetectorRadius * math.cos( deltaPhi/2.0 )
  else:
    return DetectorRadius * math.cos( deltaPhi/2.0 )


def TwoHitEvent( Event, DetectorRadius, ZMin=0.0, ZMax=0.0, RMax=120.0 ):

  if len( Event ) != 2:
    return False

  # If there's a z-cut, apply it
  if ZMin != ZMax:
    meanZ = ( Event[0][DATASET_Z] + Event[1][DATASET_Z] ) / 2.0
    if meanZ < ZMin or meanZ > ZMax:
      return False

  # Cut on the radius of closest approach
  rMin = FindHitRadius( Event, DetectorRadius )
  return math.fabs( rMin ) <= RMax

# Outdated approach
# Note that this definition specifically applies to central, linear phantoms only
def BackToBackEvent( Event, DetectorRadius, ZMin=0.0, ZMax=0.0 ):
  return TwoHitEvent( Event, DetectorRadius, ZMin, ZMax, RMax=20.0 )


# Launch the Geant4 simulation
def GenerateSample( DetectorLengthMM, Detector, SourceLengthMM, Source, TotalDecays, DetectorMaterial, Seed=1234, Path="", SourceOffset=0, NAluminiumSleeves=0 ):

  # Allow creation at arbitrary path
  outputFileName = ""
  if Path != "":
    if os.path.isdir( Path ):
      outputFileName = Path + "/"
    else:
      print( "Skipping dataset generation: unable to access path " + Path )
      return None

  # Allow for other detector granularity, but default to block
  if not ( ("Block" in Detector) or ("Panel" in Detector) or ("Crystal" in Detector) ):
      Detector += "Block"

  # Phantom length affects the attenuating material, so include it even if source is detector
  outputFileName += "hits.n" + str(TotalDecays) + "." + Detector + "." + str(DetectorLengthMM) + "mm."
  # Hanna: printing detector name commented out to have common naming convention
  # if DetectorMaterial != "":
  #   outputFileName += DetectorMaterial + "."
  outputFileName += Source + "." + str(SourceLengthMM) + "mm."
  if NAluminiumSleeves > 0 :
    outputFileName += str(NAluminiumSleeves) + "AlSleeves."
  if 'Linear' in Source:
    outputFileName += "-y" + str(SourceOffset) + "mm."

  outputFileName += str(Seed) + ".csv"

  # Check if file already present (in which case assume it's re-usable)
  if os.path.exists( outputFileName ):
    print( "Re-using previous simulation" )
    return outputFileName
  else:
    print( "Creating dataset " + outputFileName )
    command =  "../build/SimplePetScanner"
    command += " -n " + str(TotalDecays)
    command += " --detector " + Detector
    command += " --detectorLengthMM " + str(DetectorLengthMM)
    command += " --source " + Source
    command += " --phantomLengthMM " + str(SourceLengthMM)
    if NAluminiumSleeves > 0 :
      command += " --nAluminiumSleeves " + str(NAluminiumSleeves)
    command += " --outputFileName " + outputFileName
    command += " --randomSeed " + str(Seed)
    command += " --sourceOffsetMM " + str(SourceOffset)
    if DetectorMaterial != "":
      command += " --detectorMaterial " + DetectorMaterial
    print("Running command = ", command)
    process = subprocess.Popen( command, shell=True )
    process.wait()

    if process.returncode == 0:
      print( "Simulation complete" )
      return outputFileName
    else:
      print( "Simulation failed with return code: ", process.returncode )
      return ""


# Create a dataset class from new or existing simulated input
def CreateDataset( DetectorLengthMM, Detector, SourceLengthMM, Source, TotalDecays, EnergyMin, EnergyMax, DetectorMaterial, Seed=1234, Path="", ClusterLimitMM=None, SourceOffset=0, NAluminiumSleeves=0 ):

  outputFileName = GenerateSample( DetectorLengthMM, Detector, SourceLengthMM, Source, TotalDecays, DetectorMaterial, Seed, Path, SourceOffset, NAluminiumSleeves )
  if outputFileName == "":
      return None

  CLUS_DEFAULT=20
  if "Crystal" in Detector:
    if ClusterLimitMM is None:
      print( "Using a high-granularity \"Crystal\" detector geometry with no clusterisation" )
      print( "Setting to recommended value (" + str(CLUS_DEFAULT) + "mm), specify 0 to override" )
      ClusterLimitMM = CLUS_DEFAULT
    else:
      print( "Using a high-granularity \"Crystal\" detector geometry with clusterisation at " + str(ClusterLimitMM) + "mm" )

  return SimulationDataset( outputFileName, TotalDecays, EnergyMin, EnergyMax, ClusterLimitMM )
