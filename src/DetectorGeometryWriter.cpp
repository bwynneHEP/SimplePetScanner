#include "DetectorGeometryWriter.h"

#include "G4SystemOfUnits.hh"

#include <fstream>

void DetectorGeometryWriter::WriteSTIRheader( std::string const& fileName, DetectorGeometryData const& inputData )
{
  std::ofstream outputFile = std::ofstream( fileName.c_str() );

  outputFile << "ROOT header := " << std::endl;
  outputFile << std::endl;
  outputFile << "originating system := User_defined_scanner" << std::endl;
  outputFile << "Number of rings                          := " << inputData.nRings * inputData.crystalsAxial << std::endl; // Apparently each axial crystal is a ring
  outputFile << "Number of detectors per ring             := " << inputData.crystalsTrans * inputData.blocksPerRing << std::endl; // Removed axial crystals
  outputFile << "Inner ring diameter (cm)                 := " << inputData.ringInnerDiameter / cm << std::endl;
  outputFile << "Average depth of interaction (cm)        := 0.25" << std::endl; // Could calculate this value?
  outputFile << "Distance between rings (cm)              := " << inputData.ringGap / cm << std::endl; // Basically just distance between adjacent crystal centres in z I think
  outputFile << "Default bin size (cm)                    := " << inputData.crystalTransSize / cm << std::endl;
  outputFile << "View offset (degrees)                    := 0.0" << std::endl; // See STIR-UsersGuide.pdf, it's a fiddle-factor between GATE and STIR
  outputFile << "Maximum number of non-arc-corrected bins := 381" << std::endl; // What does this mean?
  outputFile << "Default number of arc-corrected bins     := 331" << std::endl; // What does this mean?
  outputFile << ";Number of TOF time bins :=275" << std::endl; // Undefined, comment out (causes STIR issue also)
  outputFile << ";Size of timing bin (ps) :=17.8" << std::endl; // ditto
  outputFile << ";Timing resolution (ps) :=75" << std::endl; // ditto
  outputFile << std::endl;
  outputFile << "GATE scanner type := GATE_Cylindrical_PET" << std::endl;
  outputFile << "GATE_Cylindrical_PET Parameters :=" << std::endl;
  outputFile << std::endl;
  outputFile << "name of data file := YOUR_FILE.root" << std::endl;
  outputFile << std::endl;
  outputFile << "name of input TChain := Coincidences" << std::endl;
  outputFile << std::endl;
  outputFile << "; As the GATE repeaters. " << std::endl;
  outputFile << "; If you skip a level in GATE's hierarchy, " << std::endl;
  outputFile << "; use 1." << std::endl;
  outputFile << "number of Rsectors := " << inputData.blocksPerRing << std::endl;
  outputFile << "number of modules_X := 1 " << std::endl;
  outputFile << "number of modules_Y := 1" << std::endl;
  outputFile << "number of modules_Z := 1" << std::endl;
  outputFile << "number of submodules_X := 1" << std::endl;
  outputFile << "number of submodules_Y := 1" << std::endl;
  outputFile << "number of submodules_Z := " << inputData.nRings << std::endl;
  outputFile << "number of crystals_X := 1" << std::endl;
  outputFile << "number of crystals_Y := " << inputData.crystalsTrans << std::endl;
  outputFile << "number of crystals_Z := " << inputData.crystalsAxial << std::endl;
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
