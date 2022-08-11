#ifndef EnergyCounter_h
#define EnergyCounter_h 1

#include "DecayTimeFinderAction.h"

#include "G4VSensitiveDetector.hh"

#include <map>
#include <fstream>

class EnergyCounter : public G4VSensitiveDetector
{
  public:
    EnergyCounter( const G4String& name, DecayTimeFinderAction * decayTimeFinder );
    ~EnergyCounter() override;

    void Initialize( G4HCofThisEvent* hitCollection ) override;
    G4bool ProcessHits( G4Step* step, G4TouchableHistory* history ) override;
    void EndOfEvent( G4HCofThisEvent* hitCollection ) override;

  private:
    std::map< G4int, G4double > m_totalEnergyMap;
    std::map< G4int, G4double > m_averageTimeMap;
    std::map< G4int, G4double > m_averageRMap;
    std::map< G4int, G4double > m_averagePhiMap;
    std::map< G4int, G4double > m_averageZMap;
    DecayTimeFinderAction * m_decayTimeFinder;
    std::ofstream m_outputFile;
};

#endif
