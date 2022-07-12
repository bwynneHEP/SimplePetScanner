#include "CrystalIntrinsicAction.h"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"

CrystalIntrinsicAction::CrystalIntrinsicAction( G4double minZ, G4double maxZ, G4double minR, G4double maxR )
  : G4VUserPrimaryGeneratorAction()
  , m_minZ( minZ ), m_maxZ( maxZ )
  , m_minR( minR ), m_maxR( maxR )
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

CrystalIntrinsicAction::~CrystalIntrinsicAction()
{
  delete m_particleGun;
}

// This function is called at the begining of event
void CrystalIntrinsicAction::GeneratePrimaries( G4Event* anEvent )
{
  // Lutetium 176
  G4int Z = 71, A = 176;
  G4double ionCharge = 0.0 * eplus;
  G4double excitEnergy = 0.0 * keV;
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon( Z, A, excitEnergy );

  // Ion at rest
  m_particleGun->SetParticleDefinition( ion );
  m_particleGun->SetParticleCharge( ionCharge );
  m_particleGun->SetParticleEnergy( 1.0 * eV );

  // Choose a position on the z-axis and a radius
  G4double z = m_minZ + ( G4UniformRand() * ( m_maxZ - m_minZ ) );
  G4double r = m_minR + ( G4UniformRand() * ( m_maxR - m_minR ) );
  G4ThreeVector position;
  position.setRhoPhiZ( r, 0.0, z );
  m_particleGun->SetParticlePosition( position );

  // Fire particle
  m_particleGun->GeneratePrimaryVertex( anEvent );
}
