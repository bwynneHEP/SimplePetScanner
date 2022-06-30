#ifndef LinearSourceAction_h
#define LinearSourceAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

// Generate a single particle and fire it into our experiment
class LinearSourceAction : public G4VUserPrimaryGeneratorAction
{
  public:
    LinearSourceAction( G4double minZ, G4double maxZ );
    ~LinearSourceAction() override;

    void GeneratePrimaries( G4Event* ) override;

  private:
    G4ParticleGun* m_particleGun;
    G4double m_minZ = 0.0;
    G4double m_maxZ = 0.0;
};

#endif
