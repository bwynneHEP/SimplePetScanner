# SimplePetScanner

This project provides basic Geant4 simulation for a PET scanner, and a set of analysis examples that use the output of the simulation.


## Installation

Requires Geant4 installation: https://geant4.web.cern.ch/

Once Geant4 is available on your system, be sure to set up the environment as follows:
```
source $GEANT4_myPath/share/Geant4-*/geant4make/geant4make.sh
```

Now clone this repository and compile (requires CMake):
```
git clone https://github.com/bwynneHEP/SimplePetScanner
cd SimplePetScanner
make
```


## Simulation

Run the compiled executable with command line arguments to set the details of the simulation, e.g.:
```
build/SimplePetScanner -n 1 --gui --detector SiemensBlock --source LinearF18
```

### -n
Set the number of events to simulate.

### --gui
Activate the gui (omit this argument to produce output data only).

### --detector
Choose a detector geometry (Siemens or Explorer) and the granularity of the detector crystals.
Panel, Block and Crystal modes are supported, corresponding to large detectors running the full axial length of the system, axial sengmentation into smaller detector modules, and fully granualar simulation at the level of individual scintillator crystals.

### --detectorLengthMM
Set the length of the detector in mm.
The actual length in simulation will be rounded up to accommodate the discrete detector modules.
Each detector geometry will take their default length unless this argument is specified.

### --detectorMaterial
Set the scintillator material used for the detector crystals.
LSO, LYSO (with the precise LYSO composition corresponding to the Explorer detector), NaI, BGO and CsF are supported.
Each detector geometry will use their default material unless this argument is specified.

### --source
Set the radioactive source, either a radioisotope tracer in a capillary (Linear) or the intrinsic detector background.
For a linear source, F18 and Zr89 isotopes are supported.
For the detector background, the detector geometry (Siemens or Explorer) should be specified, Lu176 decays will be assumed.

### --sourceOffsetMM
Shift the linear source in the negative direction of the y-axis. If this argument is non-zero and provided along the --nAluminiumSleeves argument, the sensitivity phantom is shifted as well. 

### --phantomLengthMM
Regardless of the radioactive source, a polyethylene cylindrical phantom will be simulated in the centre of the detector.
Set the length of the phantom in mm with this argument.
If a linear radioisotope source is selected, the length of the capillary will be equivalent to the length of the phantom.

### --outputFileName
Override the default output file name (hits.csv) to allow multiprocessing.

### --decayOutputFileName
Create output csv file containing the information about the radioactive decay. For now it includes only the values of the positron range for each event. The file is created only if this argument is specified in the running command. 

### --randomSeed
Override the default random seed.

### --nAluminiumSleeves
Create a sensitivity phantom with the corresponding number of aluminum sleeves.
The number of sleeves has to be an integer from 1 to 5. If the option is not called, the scatter phantom will be created by default. 

## Analysis

Example Jupyter notebooks can be found in the analysis directory, along with helper functions.
Commands that run the simulation from within the notebooks assume the relative path of the simulation executable, so please launch the notebooks from within the analysis directory.


## Acknowledgements

Includes code from https://github.com/katarinazatkova/PET-material-1
and https://gitlab.cern.ch/geant4/geant4/tree/master/examples/basic/B3
