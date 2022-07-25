# Try to reproduce result in https://www.nature.com/articles/s41598-018-35684-x

import pylab

crystalInput = open( "../crystalRadioactivity.log" )

energyMin = 300
energyMax = 600

allEnergies = []
coincidenceCounter = { 1:0, 2:0, 3:0, 4:0, 5:0 }
currentCoincidence = 1
currentEvent = -1

for line in crystalInput:

  splitLine = line.split(" ")

  energyKeV = float( splitLine[2] )
  allEnergies.append( energyKeV )

  if energyKeV < energyMin or energyKeV > energyMax:
    continue

  eventID = int( splitLine[0] )
  if currentEvent == eventID:
    currentCoincidence += 1
  else:
    if currentEvent > -1:
      coincidenceCounter[ currentCoincidence ] += 1
    currentEvent = eventID
    currentCoincidence = 1

# Finish the last coincidence (or lack thereof)
coincidenceCounter[ currentCoincidence ] += 1

print( "coincidences:", coincidenceCounter )

pylab.hist( allEnergies, bins=50 )

pylab.show()
