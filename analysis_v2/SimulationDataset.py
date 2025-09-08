# File format defined in EnergyCounter.cpp:
# EventID ModuleID[xN] Energy Time R Phi Z

# Only define positions for main quantities
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
    self.moduleMap = {}

  def __del__( self ):
    self.inputFile.close()

  def Open( self ):
    return self.nextLine != ""

  def ReadHit( self ):
    if not self.Open():
      return None

    splitLine = self.nextLine.split(" ")

    # There's a variable amount of potential module info
    moduleIDFields = len( splitLine ) - DATASET_PHOTON_LENGTH

    # Assemble hit info
    wholeHit = [ int( splitLine[DATASET_EVENT] ) ]
    for i in range( 1 + moduleIDFields, len( splitLine ) ):
      wholeHit.append( float( splitLine[i] ) )

    # Assemble module ID info
    moduleIDs = []
    for i in range( 1, moduleIDFields + 1 ):
      moduleIDs.append( int( splitLine[i] ) )

    # Save module ID map
    # Note that it's fine to use a tuple () as a dict key, but not a list []
    hitGlobalPos = ( wholeHit[ DATASET_R ], wholeHit[ DATASET_PHI ], wholeHit[ DATASET_Z ] )
    self.moduleMap[ hitGlobalPos ] = moduleIDs

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
        newEvent = self.currentEvent # save finished event
        self.currentEvent = [newHit] # open unfinished event
        return newEvent


class SimulationDataset:

  def __init__( self, InputPath, TotalDecays, EnergyMin=None, EnergyMax=None, ClusterLimitMM=None, RNG=None ):
    self.inputData = {}
    self.energyMin = EnergyMin
    self.energyMax = EnergyMax
    self.totalDecays = TotalDecays
    self.hitCount = 0
    self.eventHitsMax = 0
    self.RNG = RNG
    if self.RNG == None:
      # print("SimulationDataset.__init__ creating new RNG")
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
        self.AddHit( eventID, wholeHit, ClusterLimitMM, inputFile.moduleMap )
      else:
        # Input file should be ordered and thus old event complete
        if currentEvent in self.inputData:
          self.eventHitsMax = max( self.eventHitsMax, len( self.inputData[currentEvent] ) )

        # Start a new event
        self.inputData[ eventID ] = [ wholeHit ]
        currentEvent = eventID
        eventCount += 1
        self.hitCount += 1

    # Save the module map info
    self.moduleMap = inputFile.moduleMap

    print( str(eventCount) + " events loaded (" + str( self.totalDecays ) + " simulated) with average " + str( self.hitCount / self.totalDecays ) + " hits/event" )


  def AddHit( self, ExistingEventID, NewHit, ClusterLimitMM, ModuleMap ):

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

          # Update the module map for the merged hit
          # For now just a simple majority method: use the module with most energy
          mergedMOD = None
          if ( oldE >= newE ):
            oldHitGlobalPos = ( oldR, oldPhi, oldZ )
            mergedMOD = ModuleMap[ oldHitGlobalPos ]
          else:
            newHitGlobalPos = ( newR, newPhi, newZ )
            mergedMOD = ModuleMap[ newHitGlobalPos ]
          mergedHitGlobalPos = ( mergedR, mergedPhi, mergedZ )
          ModuleMap[ mergedHitGlobalPos ] = mergedMOD

          keepHit = False
          break

      if keepHit:
        self.inputData[ ExistingEventID ].append( NewHit )
        self.hitCount += 1


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
def CreateDataset( DetectorLengthMM, Detector, SourceLengthMM, Source, TotalDecays, EnergyMin, EnergyMax, DetectorMaterial, Seed=1234, Path="", ClusterLimitMM=None, SourceOffset=0, NAluminiumSleeves=0, UseNumpy=False ):

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

  inputData = SimulationDataset( outputFileName, TotalDecays, EnergyMin, EnergyMax, ClusterLimitMM )
  if UseNumpy:
    import NumpyDatasetReader
    return NumpyDatasetReader.NumpyDatasetReader( inputData )
  else:
    import LegacyDatasetReader
    return LegacyDatasetReader.LegacyDatasetReader( inputData )
