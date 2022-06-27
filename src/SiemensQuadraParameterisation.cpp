#include "SiemensQuadraParameterisation.h"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

SiemensQuadraParameterisation::SiemensQuadraParameterisation()
{
}

void SiemensQuadraParameterisation::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  // 32 rings in the detector
  G4int const crystalsPerRing = 7600;
  G4int const ring = floor( copyNo / crystalsPerRing );
  G4int const inRing = copyNo % crystalsPerRing;
  G4int const blocksPerRing = 38;
  //G4int const ring = floor( copyNo / blocksPerRing );
  //G4int const inRing = copyNo % blocksPerRing;

  // 38 detector blocks per ring
  G4int const crystalsPerBlock = 200;
  G4int const block = floor( inRing / crystalsPerBlock );
  G4int const inBlock = inRing % crystalsPerBlock;

  // mini-blocks are 5x5 crystals, arranged into 2x4 blocks (2 in the axial direction I think)
  // blocks are therefore 10x20
  G4int const crystalsBlockAxial = 10;
  //G4int const blockTrans = floor( inBlock / crystalsBlockAxial );
  //G4int const blockAxial = inBlock % crystalsBlockAxial;

  // Phi position is block within ring
  G4double const deltaPhi = 360.0 * deg / G4double( blocksPerRing );
  G4double const phi = deltaPhi * G4double( block );
  //G4double const phi = deltaPhi * G4double( inRing );

  // Z position is ring itself
  G4double const crystalWidth = 3.2 * mm;
  G4double const ringWidth = crystalWidth * crystalsBlockAxial;
  G4double z = G4double( ring ) * ( ringWidth + 1*mm ) - 53 * cm; // offset by half total to get centred (roughly)

  // Set the translation
  G4ThreeVector position;
  position.setRhoPhiZ( 41.0 * cm, phi, z ); // 82cm "Detector ring diameter"
  physVol->SetTranslation( position );

  // Set the rotation
  G4RotationMatrix * rotation = new G4RotationMatrix();
  rotation->rotateZ( -phi );
  physVol->SetRotation( rotation );
}
