#include "GATEheaderWriter.h"

#include <fstream>

void WriteGATEheader()
{
  std::ofstream outputFile = std::ofstream( "test.GATE.hroot" );

  outputFile << "ROOT header := " << std::endl;
  outputFile << std::endl;
  outputFile << "originating system := User_defined_scanner" << std::endl;
  outputFile << "Number of rings                          := 24" << std::endl;
  outputFile << "Number of detectors per ring             := 576" << std::endl;
  outputFile << "Inner ring diameter (cm)                 := 81.02" << std::endl;
  outputFile << "Average depth of interaction (cm)        := 0.25" << std::endl;
  outputFile << "Distance between rings (cm)              := 0.654" << std::endl;
  outputFile << "Default bin size (cm)                    := 0.21306" << std::endl;
  outputFile << "View offset (degrees)                    := -5.021" << std::endl;
  outputFile << "Maximum number of non-arc-corrected bins := 381" << std::endl;
  outputFile << "Default number of arc-corrected bins     := 331" << std::endl;
  outputFile << "Number of TOF time bins :=275" << std::endl;
  outputFile << "Size of timing bin (ps) :=17.8" << std::endl;
  outputFile << "Timing resolution (ps) :=75" << std::endl;
  outputFile << std::endl;
  outputFile << "GATE scanner type := GATE_Cylindrical_PET" << std::endl;
  outputFile << "GATE_Cylindrical_PET Parameters :=" << std::endl;
  outputFile << std::endl;
  outputFile << "name of data file := root_data_test1.root" << std::endl;
  outputFile << std::endl;
  outputFile << "name of input TChain := Coincidences" << std::endl;
  outputFile << std::endl;
  outputFile << "; As the GATE repeaters. " << std::endl;
  outputFile << "; If you skip a level in GATE's hierarchy, " << std::endl;
  outputFile << "; use 1." << std::endl;
  outputFile << "number of Rsectors := 32" << std::endl;
  outputFile << "number of modules_X := 1 " << std::endl;
  outputFile << "number of modules_Y := 1" << std::endl;
  outputFile << "number of modules_Z := 1" << std::endl;
  outputFile << "number of submodules_X := 1" << std::endl;
  outputFile << "number of submodules_Y := 2" << std::endl;
  outputFile << "number of submodules_Z := 4" << std::endl;
  outputFile << "number of crystals_X := 1" << std::endl;
  outputFile << "number of crystals_Y := 9" << std::endl;
  outputFile << "number of crystals_Z := 6" << std::endl;
  outputFile << std::endl;
  outputFile << ";; From GATE's online documentation: " << std::endl;
  outputFile << ";; (http://wiki.opengatecollaboration.org/index.php/Users_Guide_V7.2:Digitizer_and_readout_parameters)" << std::endl;
  outputFile << ";; [...] the readout depth depends upon how the electronic readout functions." << std::endl;
  outputFile << ";; If one PMT reads the four modules in the axial direction, " << std::endl;
  outputFile << ";; the depth should be set with the command:" << std::endl;
  outputFile << ";; /gate/digitizer/Singles/readout/setDepth 1 " << std::endl;
  outputFile << ";" << std::endl;
  outputFile << "; In STIR terminology this will be used to define the number of crystals" << std::endl;
  outputFile << "; per singles unit. " << std::endl;
  outputFile << "Singles readout depth := 1" << std::endl;
  outputFile << std::endl;
  outputFile << ";" << std::endl;
  outputFile << "; If set the scattered events will be skipped" << std::endl;
  outputFile << "exclude scattered events := 1" << std::endl;
  outputFile << std::endl;
  outputFile << ";" << std::endl;
  outputFile << "; If set the random events will be skipped" << std::endl;
  outputFile << "exclude random events := 1" << std::endl;
  outputFile << std::endl;
  outputFile << std::endl;
  outputFile << "; If want to deactivate set to [0, 10000]" << std::endl;
  outputFile << "low energy window (keV) := 0" << std::endl;
  outputFile << "upper energy window (keV):= 10000" << std::endl;
  outputFile << std::endl;
  outputFile << "End GATE_Cylindrical_PET Parameters :=" << std::endl;
  outputFile << std::endl;
  outputFile << "end ROOT header := " << std::endl;

  outputFile.close();
}
