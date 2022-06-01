#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4GlobalMagFieldMessenger.hh"

// Define the experiment to be simulated
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

  private:
    // Global magnetic field messenger
    static G4ThreadLocal G4GlobalMagFieldMessenger* m_magneticFieldMessenger;
};

#endif
