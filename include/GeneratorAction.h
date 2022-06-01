#ifndef GeneratorAction_h
#define GeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

// Generate a single particle and fire it into our experiment
class GeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    GeneratorAction();
    ~GeneratorAction() override;

    void GeneratePrimaries( G4Event* ) override;

  private:
    G4ParticleGun* m_particleGun;
};

#endif
