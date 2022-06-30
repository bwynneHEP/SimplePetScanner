#include "DecayTimeFinderAction.h"

#include "G4Track.hh"
//#include "G4SystemOfUnits.hh"

DecayTimeFinderAction::DecayTimeFinderAction()
{
}

DecayTimeFinderAction::~DecayTimeFinderAction()
{
}

G4ClassificationOfNewTrack DecayTimeFinderAction::ClassifyNewTrack( const G4Track* track )
{
  //std::cout << track->GetParticleDefinition()->GetParticleName() << " " << track->GetGlobalTime() - m_firstDecay << std::endl;

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
  }

  // Kill anything that took too long
  //if ( track->GetGlobalTime() - m_firstDecay > 100.0 * s ) return fKill;
  
  return fUrgent;
}
