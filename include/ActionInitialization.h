#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "DecayTimeFinderAction.h"

#include "G4VUserActionInitialization.hh"

// This class is a very simple template that just tells Geant where to look for our custom code
class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization( DecayTimeFinderAction* decayTimeFinder, std::string sourceName, G4double detectorLength, G4double phantomLength , std::string detectorMaterial, G4double sourceOffsetMM);
    ~ActionInitialization() override;

    void Build() const override;

  private:
    DecayTimeFinderAction * m_decayTimeFinder = nullptr;
    std::string m_sourceName;
    G4double m_detectorLength;
    G4double m_phantomLength;
    std::string m_detectorMaterial;
    G4double m_sourceOffsetMM;
};

#endif

