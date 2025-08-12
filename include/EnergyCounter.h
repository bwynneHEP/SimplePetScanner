#ifndef EnergyCounter_h
#define EnergyCounter_h 1

#include "DecayTimeFinderAction.h"

#include "G4VSensitiveDetector.hh"

#include <map>
#include <vector>
#include <fstream>

class EnergyCounter : public G4VSensitiveDetector
{
  public:
    EnergyCounter( const G4String& name, DecayTimeFinderAction * decayTimeFinder, std::string outputFileName );
    ~EnergyCounter() override;

    void Initialize( G4HCofThisEvent* hitCollection ) override;
    G4bool ProcessHits( G4Step* step, G4TouchableHistory* history ) override;
    void EndOfEvent( G4HCofThisEvent* hitCollection ) override;
    G4float GetEFraction( const G4int copyNo ) const;
    void SetGeometryIDs( const std::vector< std::vector< int > > GeometryIDs ){ m_geometryIDs = GeometryIDs; };

  private:
    std::map< G4int, G4double > m_totalEnergyMap;
    std::map< G4int, G4double > m_integratedEnergyMap;
    std::map< G4int, G4double > m_averageTimeMap;
    std::map< G4int, G4double > m_averageRMap;
    std::map< G4int, G4double > m_averagePhiMap;
    std::map< G4int, G4double > m_averageZMap;
    std::vector< std::vector< int > > m_geometryIDs;
    DecayTimeFinderAction * m_decayTimeFinder;
    std::ofstream m_outputFile;
    G4double m_maxEnergyValue = 0.0;
};

#endif
