#include "Parameterisation.h"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

Parameterisation::Parameterisation( G4int DetectorsPerRing ) :
  m_detectorsPerRing( DetectorsPerRing )
{
}

void Parameterisation::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  // Phi position is copyNo within ring
  G4double deltaPhi = 360.0 * deg / G4double( m_detectorsPerRing );
  G4double phi = deltaPhi * G4double( copyNo % m_detectorsPerRing );

  // Z position is ring itself
  G4int ring = floor( copyNo / m_detectorsPerRing );
  G4double z = G4double( ring ) * 20.0 * cm; // hardcoded ring spacing for now

  // Set the translation
  G4ThreeVector position;
  position.setRhoPhiZ( 50.0 * cm, phi, z ); // hardcoded radius for now
  physVol->SetTranslation( position );

  // Set the rotation
  G4RotationMatrix * rotation = new G4RotationMatrix();
  rotation->rotateZ( -phi ); // this is slow for some reason
  physVol->SetRotation( rotation );
}
