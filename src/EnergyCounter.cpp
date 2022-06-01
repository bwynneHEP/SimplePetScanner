#include "EnergyCounter.h"

#include "g4csv.hh"

EnergyCounter::EnergyCounter( const G4String& name, const G4int id )
  : G4VSensitiveDetector( name ) // Run the constructor of the parent class
{
  // Set which histogram ID / ntuple column to use
  m_ID = id;
}

EnergyCounter::~EnergyCounter()
{
}

// At the start of the event, zero the energy counter
void EnergyCounter::Initialize( G4HCofThisEvent* )
{
  m_totalEnergy = 0.0;
}

// Analyse anything that hits the detector
G4bool EnergyCounter::ProcessHits( G4Step* step, G4TouchableHistory* history )
{
  // Get the energy deposited by this hit
  G4double edep = step->GetTotalEnergyDeposit();

  // Add to the total energy in this object
  m_totalEnergy += edep;
  if ( history )
  G4cout << "BEN BEN BEN " << history->GetReplicaNumber() << G4endl;
  else
  {

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    G4cout << "BEN BEN BEN no history so use step " << theTouchable->GetReplicaNumber() << G4endl;

  }

  return true;
}

// At the end of an event, store the energy collected in this detector
void EnergyCounter::EndOfEvent( G4HCofThisEvent* )
{
  // Display the total
  G4cout << this->GetName() << " total energy = " << m_totalEnergy << G4endl;

  // Get the analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Fill histogram (histogram 0, bin by layer ID)
  analysisManager->FillH1( 0, m_ID, m_totalEnergy );

  // Fill ntuple (ntuple 0, column by layer ID)
  analysisManager->FillNtupleDColumn( 0, m_ID, m_totalEnergy );
}
