# https://github.com/UCL/STIR-GATE-Connection/blob/master/VoxelisedSimulation/SubScripts/UnlistRoot.sh
#/home/ben/Software/STIR-GATE-Connection/VoxelisedSimulation/SubScripts/UnlistRoot.sh
#/home/ben/Software/STIR-GATE-Connection/VoxelisedSimulation/DebuggingScripts/DebugUnlistRoot.sh

# Basically this means turning coincidences (list mode) into sinograms. After that I hope it's just usual STIR/SIRF stuff?
# lm_to_projdata seems to be the relevant thing, plus a par file:
#/home/ben/Software/STIR-GATE-Connection/VoxelisedSimulation/UnlistingTemplates/lm_to_projdata_template.par
# the par file relies on a .hroot header
# then just do
#lm_to_projdata THE_PAR_FILE

# variable replacements
#sed -i.bak "s|{ROOT_FILENAME}|$StoreRootFilesDirectory/${ROOT_FILENAME}|g" ${LM_TO_PROJDATA_PAR_PATH}
#sed -i.bak "s/{SinogramID}/${SinogramID}/g" ${LM_TO_PROJDATA_PAR_PATH}
#sed -i.bak "s|{UNLISTINGDIRECTORY}|${UnlistingDirectory}|g" ${LM_TO_PROJDATA_PAR_PATH}
#sed -i.bak "s|{seed}|${seed}|g" ${LM_TO_PROJDATA_PAR_PATH}
#sed -i.bak "s|{NumEventsToStore}|${NumEventsToStore}|g" ${LM_TO_PROJDATA_PAR_PATH}

# the par file stupidly also loads this template
#/home/ben/Software/STIR-GATE-Connection/ExampleScanners/mMR/STIR_scanner.hs
# which doubly-defines a lot of crap

def hrootLineParse( line, hrootDictionary ):
  if ":=" in line:
    key = line.split(":=")[0].strip()
    value = line.split(":=")[1].strip()
    hrootDictionary[key] = value

def hsLineSwap( line, hrootDictionary ):

  # Specific tweak
  if "!matrix size [3]" in line:
    numberOfViews = int( hrootDictionary["Number of detectors per ring"] ) / 2
    return "!matrix size [3] := " + str( int( numberOfViews ) ) + "\n"

  justVals = ["Number of rings", "Number of detectors per ring", "Inner ring diameter (cm)", "Average depth of interaction (cm)",
              "Distance between rings (cm)", "Default bin size (cm)", "View offset (degrees)",
              "Maximum number of non-arc-corrected bins", "Default number of arc-corrected bins" ]

  # Changes to make
  key = ""
  if "Number of crystals per block in axial direction" in line:
    key = "number of crystals_Z"
  elif "Number of crystals per block in transaxial direction" in line:
    key = "number of crystals_Y"
  elif "Number of crystals per singles unit in axial direction" in line:
    key = "Number of rings"
  elif "Number of crystals per singles unit in transaxial direction" in line:
    key = "number of crystals_Y"
  elif "Number of blocks per bucket in axial direction" in line:
    key = "number of submodules_Z"
  elif "Number of blocks per bucket in transaxial direction" in line:
    key = "number of submodules_Y"
  else:
    # Just update values
    for test in justVals:
      if test in line:
        key = test

  # Make a change if there is one
  if key == "":
    return line
  else:
    print( key )
    return line.split( ":=" )[0] + ":= " + hrootDictionary[key] + "\n"

dataFile = "/home/ben/Software/PetScanProject/SimplePetScanner/analysis_v2/testGATEoutput.root"
#hrootFile = "/home/ben/Software/PetScanProject/SimplePetScanner/TestGeomIDs/test.GATE.hroot"
hrootFile = "/home/ben/Software/PetScanProject/SimplePetScanner/TestGeomIDs/test.GATE.noTOF.hroot"
#hrootFile = "/home/ben/Software/PetScanProject/SimplePetScanner/analysis_v2/test.GATE.hroot"
newHrootFile = "test.GATE.hroot"
parFileTemplate = "/home/ben/Software/STIR-GATE-Connection/VoxelisedSimulation/UnlistingTemplates/lm_to_projdata_template.par"
newParFile = "./lm_to_projdata_template.par"
#hsFileTemplate = "/home/ben/Software/STIR-GATE-Connection/ExampleScanners/mMR/STIR_scanner.hs"
hsFileTemplate = "/home/ben/Software/PetScanProject/SimplePetScanner/Reconstruction/STIR_scanner.hs"
newHsFile = "STIR_scanner.hs"
#SinogramID = "Sino_${ROOT_FILENAME}_S${ScatterFlag}R${RandomFlag}"
sinogramID = "testSino"
unlistingDirectory = "."#testUnlist"
randomSeed = "1" # why is there a random seed?
eventsToStore = "-1" # store all

# Make a working directory
import os
topDir = os.getcwd()
jobCount = 0
while os.path.isdir( os.path.join( topDir, "job" + str(jobCount) ) ):
  jobCount += 1
jobDir = os.path.join( topDir, "job" + str(jobCount) )
os.mkdir( jobDir )
os.chdir( jobDir )

#if not os.path.isdir( unlistingDirectory ):
#  os.mkdir( unlistingDirectory )

# At the least, need to read in the hroot file, probably need to edit
hrootDictionary = {}
with open( newHrootFile, "w" ) as outputFile:
  for line in open( hrootFile ):
    if "name of data file :=" in line:
      line = "name of data file := " + dataFile
    hrootLineParse( line, hrootDictionary )
    outputFile.write( line )

# Probably need to process the hs file (to match hroot? what a shambles)
with open( newHsFile, "w" ) as outputFile:
  for line in open( hsFileTemplate ):
    line = hsLineSwap( line, hrootDictionary )
    outputFile.write( line )

# hroot gets appended automatically in STIR
if newHrootFile[-6:] == ".hroot":
  newHrootFile = newHrootFile[:-6]

# Need to process the .par file
with open( newParFile, "w" ) as outputFile:
  for line in open( parFileTemplate ):
    line = line.replace( "{ROOT_FILENAME}", newHrootFile )
    line = line.replace( "{SinogramID}", sinogramID )
    line = line.replace( "{UNLISTINGDIRECTORY}", unlistingDirectory )
    line = line.replace( "{seed}", randomSeed )
    line = line.replace( "{NumEventsToStore}", eventsToStore )
    line = line.replace( "UnlistingTemplates/STIR_scanner.hs", newHsFile )
    outputFile.write( line )

commandString = "/home/ben/Software/STIR/benInstall/bin/lm_to_projdata " + newParFile
os.system( commandString )


# seems to be the next steps
# probably need delayed too?
# Can't we just view the bloody sinogram?
#https://github.com/UCL/STIR-GATE-Connection/blob/master/DataCorrectionsComputation/ComputePoissonDataCorrections.sh
#https://github.com/UCL/STIR-GATE-Connection/blob/master/ExampleReconstruction/ExampleReconstruction.sh
