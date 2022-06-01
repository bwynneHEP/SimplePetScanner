#include "EnergyCounter.h"

#include "g4csv.hh"

EnergyCounter::EnergyCounter( const G4String& name )
  : G4VSensitiveDetector( name ) // Run the constructor of the parent class
{
}

EnergyCounter::~EnergyCounter()
{
}

// At the start of the event, zero the energy counter
void EnergyCounter::Initialize( G4HCofThisEvent* )
{
  m_totalEnergyMap.clear();
}

// Analyse anything that hits the detector
G4bool EnergyCounter::ProcessHits( G4Step* step, G4TouchableHistory* history )
{
  // Find the ID of the crystal struck this time
  G4int getID;
  if ( history )
  {
    getID = history->GetReplicaNumber();
  }
  else
  {
    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    const G4TouchableHandle& theTouchable = preStepPoint->GetTouchableHandle();
    getID = theTouchable->GetReplicaNumber();
  }

  // Get the energy deposited by this hit
  G4double edep = step->GetTotalEnergyDeposit();

  // Add to the total energy in this crystal
  m_totalEnergyMap[ getID ] += edep;

  return true;
}

// At the end of an event, store the energy collected in this detector
void EnergyCounter::EndOfEvent( G4HCofThisEvent* )
{
  // Get the analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Fill ntuple 0, column by layer ID (+1 to allow truth column)
  for ( const auto& entry : m_totalEnergyMap )
  {
    analysisManager->FillNtupleDColumn( 0, entry.first+1, entry.second );
  }
}
