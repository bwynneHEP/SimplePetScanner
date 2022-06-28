#include "EnergyCounter.h"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

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
  m_averageTimeMap.clear();
  m_averageRMap.clear();
  m_averagePhiMap.clear();
  m_averageZMap.clear();
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
  if ( edep > 0.0 )
  {
    m_totalEnergyMap[ getID ] += edep;

    // Average coordinates for energy deposit, weighted by its size
    m_averageTimeMap[ getID ] += step->GetPostStepPoint()->GetGlobalTime() * edep;
    m_averageRMap[ getID ] += step->GetPostStepPoint()->GetPosition().getRho() * edep;
    m_averagePhiMap[ getID ] += step->GetPostStepPoint()->GetPosition().getPhi() * edep;
    m_averageZMap[ getID ] += step->GetPostStepPoint()->GetPosition().z() * edep;
  }

  return true;
}

// At the end of an event, store the energy collected in this detector
void EnergyCounter::EndOfEvent( G4HCofThisEvent* )
{
  // Only output information for hits (since detector occupancy low)
  for ( const auto& entry : m_totalEnergyMap )
  {
    // Divide by the unit when outputting
    // see http://geant4.web.cern.ch/sites/geant4.web.cern.ch/files/geant4/collaboration/working_groups/electromagnetic/gallery/units/SystemOfUnits.html
    std::cout << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " " << entry.first << " " << entry.second / MeV << " ";
    std::cout << m_averageTimeMap[ entry.first ] / ( entry.second * ns ) << " ";
    std::cout << m_averageRMap[ entry.first ] / ( entry.second * mm ) << " ";
    std::cout << m_averagePhiMap[ entry.first ] / entry.second << " ";
    std::cout << m_averageZMap[ entry.first ] / ( entry.second * mm ) << std::endl;
  }
}
