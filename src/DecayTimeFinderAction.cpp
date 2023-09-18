#include "DecayTimeFinderAction.h"

#include "G4Track.hh"
#include "G4VUserTrackInformation.hh"
//#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

DecayTimeFinderAction::DecayTimeFinderAction()
{
}

DecayTimeFinderAction::~DecayTimeFinderAction()
{
}

G4ClassificationOfNewTrack DecayTimeFinderAction::ClassifyNewTrack( const G4Track* track )
{
  

  // Reset the timer if the track has no parent (i.e. it's a new event)
  if ( track->GetParentID() == 0 )
  {
    m_firstDecay = 0.0;
    //std::cout << "RESET" << std::endl;
  }

  // Find the first track with non-zero time
  if ( m_firstDecay == 0.0 && track->GetGlobalTime() != 0.0 )
  {
    m_firstDecay = track->GetGlobalTime();
    
    //std::cout << "SET DECAY TIME" << std::endl;

    //Get the position of the radioactive decay
    m_radDecayX = track->GetPosition().x();
    m_radDecayY = track->GetPosition().y();
    m_radDecayZ = track->GetPosition().z();
  }

  // Get the position of the electron-positron annihilation
  if (track->GetParticleDefinition()->GetParticleName() == "gamma" && track->GetParentID() == 4){
    G4ThreeVector annihilationPos = track->GetPosition();
    m_annihilationX = annihilationPos.x();
    m_annihilationY = annihilationPos.y();
    m_annihilationZ = annihilationPos.z();
  }
  // Kill anything that took too long
  //if ( track->GetGlobalTime() - m_firstDecay > 100.0 * s ) return fKill;
  
  return fUrgent;
}

G4double DecayTimeFinderAction::GetPositronRange()
{
  G4double d2 = pow(m_annihilationX-m_radDecayX, 2) + pow(m_annihilationY-m_radDecayY, 2) + pow(m_annihilationZ-m_radDecayZ, 2);
  m_positronRange = sqrt(d2);
  return m_positronRange;
}