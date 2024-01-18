#ifndef DecayTimeFinderAction_h
#define DecayTimeFinderAction_h 1

#include "G4UserStackingAction.hh"
#include "G4Types.hh"

#include <fstream>
#include <vector>

// Finds the first new track with a non-zero global time
class DecayTimeFinderAction : public G4UserStackingAction
{
  public:
    DecayTimeFinderAction(std::string decayOutputFileName);
    ~DecayTimeFinderAction() override;
     
    G4ClassificationOfNewTrack ClassifyNewTrack( const G4Track* ) override;

    void NewStage() override;

    G4double GetDecayTime() const
    {
      return m_firstDecay;
    }

    std::vector<G4double> GetDecayPos()
    {
      return {m_radDecayX, m_radDecayY, m_radDecayZ};
    }

    std::vector<G4double> GetAnnihilationPos()
    {
      return {m_annihilationX, m_annihilationY, m_annihilationZ};
    }

    G4double GetPositronRange();

  private:
    G4double m_firstDecay = 0.0;
    G4double m_radDecayX = 0.0;
    G4double m_radDecayY = 0.0;
    G4double m_radDecayZ = 0.0;
    G4double m_annihilationX = 0.0;
    G4double m_annihilationY = 0.0;
    G4double m_annihilationZ = 0.0;
    std::ofstream m_decayOutputFile;
};

#endif

