import os

def HROOTlineParse( line, HROOTdictionary ):
  if ":=" in line:
    key = line.split(":=")[0].strip()
    value = line.split(":=")[1].strip()
    HROOTdictionary[key] = value

def HSlineSwap( line, HROOTdictionary ):

  # Specific calculation needs doing
  if "!matrix size [3]" in line:
    numberOfViews = int( HROOTdictionary["Number of detectors per ring"] ) / 2
    return "!matrix size [3] := " + str( int( numberOfViews ) ) + "\n"
  
  # HACK - try to understand how to define correctly
  if "!matrix size [4]" in line:
    return "!matrix size [4] := 1\n"
  if "!matrix size [2]" in line:
    #return "!matrix size [2] := { 1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001 }\n"
    #2001 -> 4954170 stored, in sinograms 1436:1847 # increase ring difference and events up: 12423700 (in same sinos)
    #1001 -> 645519 stored, in sinograms 936:1000
    #501 -> 0 stored
    return "!matrix size [2] := { 2001 }\n"
  if "minimum ring difference per segment" in line:
    return "minimum ring difference per segment := { -320 }\n" #-60
  if "maximum ring difference per segment" in line:
    return "maximum ring difference per segment := { 320 }\n" #60

  # Things to copy directly
  justVals = ["Number of rings", "Number of detectors per ring", "Inner ring diameter (cm)", "Average depth of interaction (cm)",
              "Distance between rings (cm)", "Default bin size (cm)", "View offset (degrees)",
              "Maximum number of non-arc-corrected bins", "Default number of arc-corrected bins" ]

  # Mapping of equivalent values with different names
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
    # Just copy values with same names
    for test in justVals:
      if test in line:
        key = test

  # Make a change if there is one requested
  if key == "":
    return line
  else:
    print( key )
    return line.split( ":=" )[0] + ":= " + HROOTdictionary[key] + "\n"


def CreateSinogram( InputFileName, OutputSinogramName, HSfileTemplate="STIR_scanner.hs", PARfileTemplate="lm_to_projdata_template.par" ):

  # File names
  HROOTfileName = InputFileName + ".hroot"
  ROOTfileName = InputFileName + ".root"
  HSfileName = InputFileName + ".hs"
  PARfileName = InputFileName + ".par"

  # Parse the HROOT file
  # Overwrite to point to the correct ROOT file
  HROOTdictionary = {}
  HROOTtext = ""
  with open( HROOTfileName, "r" ) as inputHROOTfile:
    HROOTtext = inputHROOTfile.read() # Just store the file in memory
  with open( HROOTfileName, "w" ) as outputHROOTfile:
    for line in HROOTtext.split( "\n" ):
      if "name of data file :=" in line:
        line = "name of data file := " + ROOTfileName
      HROOTlineParse( line, HROOTdictionary )
      outputHROOTfile.write( line + "\n" ) # Overwrite the file

  # Create an HS file with duplicate information from HROOT
  with open( HSfileName, "w" ) as outputHSfile:
    for line in open( HSfileTemplate, "r" ):
      line = HSlineSwap( line, HROOTdictionary )
      outputHSfile.write( line )

  # Create a PAR file that just stores task parameters
  with open( PARfileName, "w" ) as outputPARfile:
    for line in open( PARfileTemplate, "r" ):
      line = line.replace( "{ROOT_FILENAME}", InputFileName ) # .hroot is appended internally
      line = line.replace( "{SinogramID}", OutputSinogramName )
      line = line.replace( "{UNLISTINGDIRECTORY}", "." ) # operate in local directory
      line = line.replace( "{seed}", "1" ) # why do we need a random seed?
      line = line.replace( "{NumEventsToStore}", "-1" ) # store all events
      line = line.replace( "UnlistingTemplates/STIR_scanner.hs", HSfileName )
      outputPARfile.write( line )

  # Assumes STIR executables are in $PATH
  commandString = "lm_to_projdata " + PARfileName
  os.system( commandString )
