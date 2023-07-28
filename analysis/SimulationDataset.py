# File format defined in EnergyCounter.cpp:
#   std::cout << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " " << entry.first << " " << entry.second / keV << " ";
#   std::cout << m_averageTimeMap[ entry.first ] / ( entry.second * ns ) << " ";
#   std::cout << m_averageRMap[ entry.first ] / ( entry.second * mm ) << " ";
#   std::cout << m_averagePhiMap[ entry.first ] / entry.second << " ";
#   std::cout << m_averageZMap[ entry.first ] / ( entry.second * mm ) << std::endl;
# so
# EventID ModuleID Energy Time R Phi Z

import random
import math
import subprocess
import os
import numpy as np

class SimulationDataset:

  def __init__( self, InputPath, TotalDecays, EnergyMin, EnergyMax, CoincidenceWindow=0.0 ):
    self.inputData = {}
    self.unusedEvents = []
    self.usedEvents = []
    self.totalDecays = TotalDecays

    CoincidenceWindow *= 1E9 #nanoseconds
    eventCount = 0.0
    hitCount = 0.0
    currentEvent = -1
    eventStartTime = 0.0

    # Parse input
    inputFile = open( InputPath )
    for line in inputFile:
      splitLine = line.split(" ")

      eventID = int( splitLine[0] )
      moduleID = int( splitLine[1] )
      energyKeV = float( splitLine[2] )
      timeNS = float( splitLine[3] )

      # If the file is ordered (it should be) then events will be contiguous
      if eventID == currentEvent:
        # Find the earliest time in the current event
        eventStartTime = min( eventStartTime, timeNS )
      else:
        # Apply the time cut to the old event
        if CoincidenceWindow > 0.0 and currentEvent in self.inputData:
          for hit in self.inputData[ currentEvent ]:
            if hit[2] > eventStartTime + CoincidenceWindow:
              self.inputData[ currentEvent ].remove( hit )

        # Start a new event
        eventStartTime = timeNS
        currentEvent = eventID

      # Do the energy window on loading to keep RAM down
      if energyKeV >= EnergyMin and energyKeV <= EnergyMax:
        wholeEvent = [ eventID, moduleID, energyKeV, timeNS ]
        for i in range( 4, len( splitLine ) ):
          wholeEvent.append( float( splitLine[i] ) )

        # Multiple lines (hits) can go into a single event
        if eventID in self.inputData:
          self.inputData[ eventID ].append( wholeEvent )
        else:
          self.inputData[ eventID ] = [ wholeEvent ]
          eventCount += 1.0
        hitCount += 1.0

    inputFile.close()

    if eventCount > 0:
      print( str(eventCount) + " events loaded (" + str( self.totalDecays ) + " simulated) with average " + str( hitCount / eventCount ) + " hits/event" )
    else:
      print( str(eventCount) + " events loaded (" + str( self.totalDecays ) + " simulated)" )

    # Allow for decays that don't enter the energy window
    for i in range( self.totalDecays ):
      self.unusedEvents.append( i )

  def SampleOneEvent( self, EnergyResolution=0.0, TimeResolution=0.0 ):
    # Check if we have any events left
    if len( self.unusedEvents ) == 0:
      self.unusedEvents = self.usedEvents
      random.shuffle( self.unusedEvents )
      self.usedEvents = []

    eventID = self.unusedEvents.pop(-1)
    self.usedEvents.append( eventID )
    if eventID in self.inputData:
      if EnergyResolution > 0.0 or TimeResolution > 0.0:

        modifiedEvent = []
        for photon in self.inputData[ eventID ]:

          newEnergy = photon[1] * ( 1 + np.random.normal( 0.0, EnergyResolution ) ) # Energy resolution as a percentage
          newTime = photon[2] + ( np.random.normal( 0.0, TimeResolution ) ) # Time resolution as absolute ns
          modifiedEvent += [[ photon[0], newEnergy, newTime, photon[3], photon[4], photon[5] ]]

        return modifiedEvent

      else:
        return self.inputData[ eventID ]
    else:
      return []

  def size( self ):
    return self.totalDecays

def FindHitRadius( Event, DetectorRadius):
  if len( Event ) != 2:
    return -1.0

  #check event ID from both photons 
  hasSameEventID = False
  if (Event[0][0] == Event[1][0]) :
    hasSameEventID = True

  # Calculate delta phi
  phi1 = Event[0][5]
  phi2 = Event[1][5]
  deltaPhi = phi1 - phi2
  while deltaPhi > math.pi:
    deltaPhi -= 2.0 * math.pi
  while deltaPhi < -math.pi:
    deltaPhi += 2.0 * math.pi

  if deltaPhi < 0.0:
    return -DetectorRadius * math.cos( deltaPhi/2.0 ), hasSameEventID
  else:
    return DetectorRadius * math.cos( deltaPhi/2.0 ), hasSameEventID

def TwoHitEvent( Event, DetectorRadius, ZMin=0.0, ZMax=0.0, RMax=120.0 ):
  if len( Event ) != 2:
    return False

  # If there's a z-cut, apply it
  if ZMin != ZMax:
    meanZ = ( Event[0][6] + Event[1][6] ) / 2.0
    if meanZ < ZMin or meanZ > ZMax:
      return False

  # Cut on the radius of closest approach
  rMin, hasSameEventID = FindHitRadius( Event, DetectorRadius )
  return math.fabs( rMin ) <= RMax

def BackToBackEvent( Event, DetectorRadius, ZMin=0.0, ZMax=0.0 ):
  return TwoHitEvent( Event, DetectorRadius, ZMin, ZMax, RMax=20.0 )

def CreateDataset( DetectorLengthMM, Detector, SourceLengthMM, Source, TotalDecays, EnergyMin, EnergyMax, DetectorMaterial, Seed=1234, CoincidenceWindow=0.0 ):
  print("DetectorMaterial in createdataset: ", DetectorMaterial)
  # Phantom length affects the attenuating material, so include it even if source is detector
  outputFileName = "hits.n" + str(TotalDecays) + "." + Detector + "Block." + str(DetectorLengthMM) + "mm."
  # Hanna: printing detector name commented out to have common naming convention
  # if DetectorMaterial != "":
  #   outputFileName += DetectorMaterial + "."
  outputFileName += Source + "." + str(SourceLengthMM) + "mm." + str(Seed) + ".csv"

  # Check if file already present (in which case assume it's re-usable)
  if os.path.exists( outputFileName ):
    print( "Re-using previous simulation" )
  else:
    command =  "../build/SimplePetScanner"
    command += " -n " + str(TotalDecays)
    command += " --detector " + Detector + "Block"
    command += " --detectorLengthMM " + str(DetectorLengthMM)
    command += " --source " + Source
    command += " --phantomLengthMM " + str(SourceLengthMM)
    command += " --outputFileName " + outputFileName
    command += " --randomSeed " + str(Seed)
    if DetectorMaterial != "":
      command += " --detectorMaterial " + DetectorMaterial
    print("Running command = ", command)
    process = subprocess.Popen( command, shell=True )
    process.wait() # Later can do some multiprocess stuff
    if process.returncode == 0:
      print( "Simulation complete" )
    else:
      print( "Simulation failed with return code: ", process.returncode )
      return None
  return SimulationDataset( outputFileName, TotalDecays, EnergyMin, EnergyMax, CoincidenceWindow )
