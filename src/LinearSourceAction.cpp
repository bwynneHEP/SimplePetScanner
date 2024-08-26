#include "LinearSourceAction.h"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"

#include <utility>

LinearSourceAction::LinearSourceAction( G4double minZ, G4double maxZ, std::string isotope )
  : G4VUserPrimaryGeneratorAction()
  , m_minZ( minZ ), m_maxZ( maxZ )
  , m_isotope( isotope )
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
  if ( m_isotope.size() )
  {
    G4int Z = 0, A = 0;

    std::map<std::string, std::pair<G4int, G4int>> isotopeMap = {
      {"F18", {9, 18}},
      {"Zr89", {40, 89}},
      {"Y90", {39, 90}},
      {"C11", {6, 11}}, 
      {"O15", {8, 15}},
      {"N13", {7, 13}},
      {"Rb82", {37, 82}},
      {"Ga68", {31, 68}},
    };

    if (auto search = isotopeMap.find(m_isotope); search != isotopeMap.end()){
      Z = search->second.first;
      A = search->second.second;
    }
    else {
      std::cerr << "Unrecognised tracer isotope: " << m_isotope << std::endl;
      exit(1);
    }

    G4double ionCharge = 0.0 * eplus;
    G4double excitEnergy = 0.0 * keV;
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon( Z, A, excitEnergy );

    // Ion at rest
    m_particleGun->SetParticleDefinition( ion );
    m_particleGun->SetParticleCharge( ionCharge );
    m_particleGun->SetParticleEnergy( 1.0 * eV );
  }

  // Choose a position on the z-axis
  G4double z = m_minZ + ( G4UniformRand() * ( m_maxZ - m_minZ - (40*mm) ) );
  //G4double z = m_minZ + ( (1.0/G4UniformRand()) * ( m_maxZ - m_minZ ) );
  G4double r = G4UniformRand() * 20*mm;
  G4double x = r * G4UniformRand();
  G4double y = sqrt( (r*r) - (x*x) );
  m_particleGun->SetParticlePosition( G4ThreeVector( x, y, z ) );

  // Fire particle
  m_particleGun->GeneratePrimaryVertex( anEvent );
}
