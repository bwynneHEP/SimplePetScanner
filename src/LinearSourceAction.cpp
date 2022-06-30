#include "LinearSourceAction.h"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"

LinearSourceAction::LinearSourceAction( G4double minZ, G4double maxZ )
  : G4VUserPrimaryGeneratorAction()
  , m_minZ( minZ ), m_maxZ( maxZ )
{
  G4int nofParticles = 1;
  m_particleGun = new G4ParticleGun( nofParticles );

  // 511 keV photon (default, probably overridden later)
  G4ParticleDefinition * particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle( "gamma" );
  m_particleGun->SetParticleDefinition( particleDefinition );
  m_particleGun->SetParticleEnergy( 511*keV );

  m_particleGun->SetParticlePosition( G4ThreeVector( 0.0, 0.0, 0.0 ) ); // in the middle of the detector
  m_particleGun->SetParticleMomentumDirection( G4ThreeVector( 0.0, 0.0, 1.0 ) ); // along z axis
}

LinearSourceAction::~LinearSourceAction()
{
  delete m_particleGun;
}

// This function is called at the begining of event
void LinearSourceAction::GeneratePrimaries( G4Event* anEvent )
{
  // Fluorine 18
  G4int Z = 9, A = 18;
  G4double ionCharge = 0.0 * eplus;
  G4double excitEnergy = 0.0 * keV;
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon( Z, A, excitEnergy );

  // Zirconium 89
/*  G4int Z = 40, A = 89;
  G4double ionCharge = 0.0 * eplus;
  G4double excitEnergy = 0.0 * keV;
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon( Z, A, excitEnergy );*/

  // Yttrium 90
/*  G4int Z = 39, A = 90;
  G4double ionCharge = 0.0 * eplus;
  G4double excitEnergy = 0.0 * keV;
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon( Z, A, excitEnergy );*/

  // Ion at rest
  m_particleGun->SetParticleDefinition( ion );
  m_particleGun->SetParticleCharge( ionCharge );
  m_particleGun->SetParticleEnergy( 1.0 * eV );

  // Choose a position on the z-axis
  G4double z = m_minZ + ( G4UniformRand() * ( m_maxZ - m_minZ ) );
  m_particleGun->SetParticlePosition( G4ThreeVector( 0.0, 0.0, z ) );

  // Fire particle
  m_particleGun->GeneratePrimaryVertex( anEvent );
}
