#include "EnergyCounter.h"

#include "G4RunManager.hh"

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
  if ( edep > 0.0 ) m_totalEnergyMap[ getID ] += edep;

  return true;
}

// At the end of an event, store the energy collected in this detector
void EnergyCounter::EndOfEvent( G4HCofThisEvent* )
{
  // Only output information for hits (since detector occupancy low)
  for ( const auto& entry : m_totalEnergyMap )
  {
    std::cout << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " " << entry.first << " " << entry.second << std::endl;
  }
}
