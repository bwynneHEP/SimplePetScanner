#include "SiemensQuadraParameterisationPanels.h"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

SiemensQuadraParameterisationPanels::SiemensQuadraParameterisationPanels()
{
}

void SiemensQuadraParameterisationPanels::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  // 32 rings in the detector
  G4int const blocksPerRing = 38;

  // Phi position is block within ring
  G4double const deltaPhi = 360.0 * deg / G4double( blocksPerRing );
  G4double const phi = deltaPhi * G4double( copyNo );

  G4double const z = 0.0; // panels go the full length

  G4double const r = 41.0 * cm; // 82cm "Detector ring diameter"

  // Set the translation
  G4ThreeVector position;
  position.setRhoPhiZ( r, phi, z );
  physVol->SetTranslation( position );

  // Set the rotation
  G4RotationMatrix * rotation = new G4RotationMatrix();
  rotation->rotateZ( -phi );
  physVol->SetRotation( rotation );
}
