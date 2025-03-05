#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "DecayTimeFinderAction.h"
#include "EnergyCounter.h"

#include "G4VUserDetectorConstruction.hh"
#include "G4GlobalMagFieldMessenger.hh"

// Define the experiment to be simulated
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction( DecayTimeFinderAction * decayTimeFinder, std::string detector, G4double detectorLength, G4double phantomLength, std::string outputFileName, std::string decayOutputFileName, std::string material, G4int nAluminiumSleeves );
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

  private:
    // Global magnetic field messenger
    static G4ThreadLocal G4GlobalMagFieldMessenger* m_magneticFieldMessenger;

    DecayTimeFinderAction * m_decayTimeFinder;
    EnergyCounter * m_energyCounter;
    std::string m_detector;
    std::string m_outputFileName;
    std::string m_decayOutputFileName;
    std::string m_material;
    G4double m_detectorLength;
    G4double m_phantomLength;
    G4int m_nAluminiumSleeves;
};

#endif
