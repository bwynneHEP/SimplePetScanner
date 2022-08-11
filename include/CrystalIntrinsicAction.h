#ifndef CrystalIntrinsicAction_h
#define CrystalIntrinsicAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

// Generate a single particle and fire it into our experiment
class CrystalIntrinsicAction : public G4VUserPrimaryGeneratorAction
{
  public:
    CrystalIntrinsicAction( G4double minZ, G4double maxZ, std::string crystalType, G4double minR, G4double maxR );
    ~CrystalIntrinsicAction() override;

    void GeneratePrimaries( G4Event* ) override;

  private:
    G4ParticleGun* m_particleGun;
    G4double m_minZ = 0.0;
    G4double m_maxZ = 0.0;
    G4double m_minR = 0.0;
    G4double m_maxR = 0.0;
    std::string m_crystalType;
};

#endif
