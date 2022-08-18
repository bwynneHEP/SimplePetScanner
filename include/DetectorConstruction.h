#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "DecayTimeFinderAction.h"

#include "G4VUserDetectorConstruction.hh"
#include "G4GlobalMagFieldMessenger.hh"

// Define the experiment to be simulated
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction( DecayTimeFinderAction * decayTimeFinder, std::string detector, G4double detectorLength, G4double phantomLength, std::string outputFileName );
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

  private:
    // Global magnetic field messenger
    static G4ThreadLocal G4GlobalMagFieldMessenger* m_magneticFieldMessenger;

    DecayTimeFinderAction * m_decayTimeFinder;
    std::string m_detector;
    std::string m_outputFileName;
    G4double m_detectorLength;
    G4double m_phantomLength;
};

#endif
