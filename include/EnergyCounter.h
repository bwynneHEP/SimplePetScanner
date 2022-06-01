#ifndef EnergyCounter_h
#define EnergyCounter_h 1

#include "G4VSensitiveDetector.hh"


class EnergyCounter : public G4VSensitiveDetector
{
  public:
    EnergyCounter( const G4String& name, const G4int id );
    ~EnergyCounter() override;

    void Initialize( G4HCofThisEvent* hitCollection ) override;
    G4bool ProcessHits( G4Step* step, G4TouchableHistory* history ) override;
    void EndOfEvent( G4HCofThisEvent* hitCollection ) override;

  private:
    G4double m_totalEnergy;
    G4int m_ID;
};

#endif
