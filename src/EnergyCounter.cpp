#include "EnergyCounter.h"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

EnergyCounter::EnergyCounter( const G4String& name, DecayTimeFinderAction * decayTimeFinder, std::string outputFileName )
  : G4VSensitiveDetector( name ) // Run the constructor of the parent class
  , m_decayTimeFinder( decayTimeFinder )
  , m_outputFile( outputFileName )
{
  if (outputFileName.empty()){
    std::cerr << "Output file name cannot be empty" << std::endl;
    exit(1);
  }
  if (!m_outputFile.good()){
    std::cerr << "Failed to open file: " << outputFileName << std::endl;
    exit(1);
  }
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
    m_integratedEnergyMap[ getID ] += edep;
    if ( m_integratedEnergyMap[ getID ] > m_maxEnergyValue && getID ) m_maxEnergyValue = m_integratedEnergyMap[ getID ];

    // Average coordinates for energy deposit, weighted by its size
    m_averageTimeMap[ getID ] += ( step->GetPostStepPoint()->GetGlobalTime() - m_decayTimeFinder->GetDecayTime() ) * edep;
    m_averageRMap[ getID ] += step->GetPostStepPoint()->GetPosition().getRho() * edep;
    m_averagePhiMap[ getID ] += step->GetPostStepPoint()->GetPosition().getPhi() * edep;
    m_averageZMap[ getID ] += step->GetPostStepPoint()->GetPosition().z() * edep;
  }

  return true;
}

// At the end of an event, store the energy collected in this detector
void EnergyCounter::EndOfEvent( G4HCofThisEvent* )
{
  // Check if the geometry ID lookup is compatible with this data
  int maxID = 0;
  for ( const auto& entry : m_totalEnergyMap )
  {
    if ( entry.first > maxID ) maxID = entry.first;
  }
  bool useGeometryIDs = ( maxID < m_geometryIDs.size() );

  // Only output information for hits (since detector occupancy low)
  for ( const auto& entry : m_totalEnergyMap )
  {
    m_outputFile << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " ";

    // If there is a geometry ID lookup (e.g. for STIR/GATE compatibility) then use it
    if ( useGeometryIDs )
    {
      // Can be multiple ID values, e.g. crystal, block, ring, etc.
      for ( const auto& idValue : m_geometryIDs[ entry.first ] ) m_outputFile << idValue << " ";
    }
    else
    {
      m_outputFile << entry.first << " ";
    }

    // Divide by the unit when outputting
    // see http://geant4.web.cern.ch/sites/geant4.web.cern.ch/files/geant4/collaboration/working_groups/electromagnetic/gallery/units/SystemOfUnits.html
    m_outputFile << entry.second / keV << " ";
    m_outputFile << m_averageTimeMap[ entry.first ] / ( entry.second * ns ) << " ";
    m_outputFile << m_averageRMap[ entry.first ] / ( entry.second * mm ) << " ";
    m_outputFile << m_averagePhiMap[ entry.first ] / entry.second << " ";
    m_outputFile << m_averageZMap[ entry.first ] / ( entry.second * mm ) << std::endl;
  }
}

G4float EnergyCounter::GetEFraction( const G4int copyNo ) const
{
  auto searchResult = m_integratedEnergyMap.find( copyNo );
  if ( searchResult != m_integratedEnergyMap.end() )
  {
    // Use a log scale, helps for the intrinsic activity
    double logVal = log( searchResult->second / m_maxEnergyValue );
    return 1.0 + ( logVal * 0.1 ); // logs will be negative
  }
  return 0.0;
}
