#ifndef DecayTimeFinderAction_h
#define DecayTimeFinderAction_h 1

#include "G4UserStackingAction.hh"
#include "G4Types.hh"

// Finds the first new track with a non-zero global time
class DecayTimeFinderAction : public G4UserStackingAction
{
  public:
    DecayTimeFinderAction();
    ~DecayTimeFinderAction() override;
     
    G4ClassificationOfNewTrack ClassifyNewTrack( const G4Track* ) override;

    G4double GetDecayTime() const
    {
      return m_firstDecay;
    }

  private:
    G4double m_firstDecay = 0.0;
};

#endif

