#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "DecayTimeFinderAction.h"

#include "G4VUserDetectorConstruction.hh"
#include "G4GlobalMagFieldMessenger.hh"

// Define the experiment to be simulated
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction( DecayTimeFinderAction * decayTimeFinder );
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

  private:
    // Global magnetic field messenger
    static G4ThreadLocal G4GlobalMagFieldMessenger* m_magneticFieldMessenger;

    DecayTimeFinderAction * m_decayTimeFinder;
};

#endif
