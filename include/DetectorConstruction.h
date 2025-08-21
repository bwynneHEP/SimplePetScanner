#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "DecayTimeFinderAction.h"
#include "EnergyCounter.h"
#include "DetectorGeometryWriter.h"

#include "G4VUserDetectorConstruction.hh"
#include "G4GlobalMagFieldMessenger.hh"

// Define the experiment to be simulated
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction( DecayTimeFinderAction * decayTimeFinder, std::string detector, G4double detectorLength, G4double phantomLength, std::string outputFileName, std::string material, G4int nAluminiumSleeves, G4double sourceOffsetMM );
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    void MakeHeaderFile( std::string const& headerFileName );

  private:
    // Global magnetic field messenger
    static G4ThreadLocal G4GlobalMagFieldMessenger* m_magneticFieldMessenger;

    DecayTimeFinderAction * m_decayTimeFinder;
    EnergyCounter * m_energyCounter;
    std::string m_detector;
    std::string m_outputFileName;
    std::string m_material;
    G4double m_detectorLength;
    G4double m_phantomLength;
    G4int m_nAluminiumSleeves;
    G4double m_sourceOffsetMM;
    DetectorGeometryData m_detectorData;
};

#endif
