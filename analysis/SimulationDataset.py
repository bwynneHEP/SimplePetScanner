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

class SimulationDataset:

  def __init__( self, InputPath, TotalDecays, EnergyMin, EnergyMax ):
    self.inputData = {}
    self.unusedEvents = []
    self.usedEvents = []

    eventCount = 0.0
    hitCount = 0.0

    # Parse input
    inputFile = open( InputPath )
    for line in inputFile:
      splitLine = line.split(" ")

      eventID = int( splitLine[0] )
      moduleID = int( splitLine[1] )
      energyKeV = float( splitLine[2] )

      # Do the energy window on loading to keep RAM down
      if energyKeV >= EnergyMin and energyKeV <= EnergyMax:
        wholeEvent = [ moduleID, energyKeV ]
        for i in range( 3, len( splitLine ) ):
          wholeEvent.append( float( splitLine[i] ) )

        # Multiple lines (hits) can go into a single event
        if eventID in self.inputData:
          self.inputData[ eventID ].append( wholeEvent )
        else:
          self.inputData[ eventID ] = [ wholeEvent ]
          eventCount += 1.0
        hitCount += 1.0

    inputFile.close()

    print( str(eventCount) + " events loaded (" + str(TotalDecays) + " simulated) with average " + str( hitCount / eventCount ) + " hits/event" )

    # Allow for decays that don't enter the energy window
    self.totalDecays = TotalDecays
    for i in range( self.totalDecays ):
      self.unusedEvents.append( i )

  def SampleOneEvent( self ):

    # Check if we have any events left
    if len( self.unusedEvents ) == 0:
      self.unusedEvents = self.usedEvents
      random.shuffle( self.unusedEvents )
      self.usedEvents = []

    eventID = self.unusedEvents.pop(-1)
    self.usedEvents.append( eventID )
    if eventID in self.inputData:
      return self.inputData[ eventID ]
    else:
      return []

def TwoHitEvent( Event ):
  return len( Event ) == 2

def BackToBackEvent( Event, PhiTolerance ):
  if not TwoHitEvent( Event ):
    return False

  phi1 = Event[0][4]
  phi2 = Event[1][4]

  deltaPhi = phi1 - phi2
  while deltaPhi > math.pi:
    deltaPhi -= 2.0 * math.pi
  while deltaPhi < -math.pi:
    deltaPhi += 2.0 * math.pi
  return math.fabs( deltaPhi ) >= ( math.pi - PhiTolerance )

def CreateDataset( LengthMM, Source, TotalDecays, EnergyMin, EnergyMax ):
  outputFileName = "hits.n" + str(TotalDecays) + ".SiemensBlock." + Source + "." + str(LengthMM) + "mm.csv"

  # Check if file already present (in which case assume it's re-usable)
  if not os.path.exists( outputFileName ):
    command =  "../build/SimplePetScanner"
    command += " -n " + str(TotalDecays)
    command += " --detector SiemensBlock"
    command += " --source " + Source
    command += " --detectorLengthMM " + str(LengthMM)
    command += "; mv hits.csv " + outputFileName
    process = subprocess.Popen( command, shell=True )
    process.wait() # Later can do some multiprocess stuff if fix the file names

  return SimulationDataset( outputFileName, TotalDecays, EnergyMin, EnergyMax )
