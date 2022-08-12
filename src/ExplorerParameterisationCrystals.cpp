#include "ExplorerParameterisationCrystals.h"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

ExplorerParameterisationCrystals::ExplorerParameterisationCrystals( G4int nCopies )
{
  // Precalculating everything avoids a memory leak
  m_positions.reserve( nCopies );
  m_rotations.reserve( nCopies );
  m_visions.reserve( nCopies );

  for ( G4int copyNo = 0; copyNo < nCopies; ++copyNo )
  {
    // 32 rings in the detector
    G4int const crystalsPerRing = 70560;
    G4int const blocksPerRing = 24;
    G4int const ring = floor( copyNo / crystalsPerRing );
    G4int const nRings = ceil( nCopies / crystalsPerRing );
    G4int const inRing = copyNo % crystalsPerRing;

    // 38 detector blocks per ring
    G4int const crystalsPerBlock = 2940;
    G4int const block = floor( inRing / crystalsPerBlock );
    G4int const inBlock = inRing % crystalsPerBlock;

    // mini-blocks are 6x7 crystals, arranged into 14x5 blocks
    // blocks are therefore 84x35
    G4int const crystalsBlockAxial = 84;
    G4int const crystalsBlockTrans = 35;
    G4int const blockTrans = floor( inBlock / crystalsBlockAxial );
    G4int const blockAxial = inBlock % crystalsBlockAxial;

    // Phi position is block within ring
    G4double const deltaPhi = 360.0 * deg / G4double( blocksPerRing );
    G4double const phi = deltaPhi * G4double( block );

    G4double const ringOffset = 2.6*mm;

    // Z position is ring itself
    G4double const crystalWidth = 2.76 * mm;
    G4double const ringWidth = crystalWidth * crystalsBlockAxial;
    //G4double const z = ( G4double( ring ) - G4double( nRings - 1 ) / 2.0 ) * ( ringWidth + 1*mm ); // +1mm for easy view
    G4double const z = ( G4double( ring ) - G4double( nRings - 1 ) / 2.0 ) * ( ringWidth + ringOffset ); // offset to get it centred

    // Adjust Z position for the block axial
    //G4double const dZ = ( crystalWidth + 1*mm ) * ( blockAxial - ( crystalsBlockAxial / 2.0 ) ); // +1mm for easy view
    G4double const dZ = crystalWidth * ( blockAxial - ( crystalsBlockAxial / 2.0 ) );

    // Adjust phi position for block transaxial
    G4double const r = 40.205 * cm; // 78.6 cm in diameter (detector face-to-face) + 18.1/2 crystal depth
    //G4double const ta = ( crystalWidth + 1*mm ) * ( blockTrans - ( crystalsBlockTrans / 2.0 ) ); // +1mm for easy view
    G4double const ta = crystalWidth * ( blockTrans - ( crystalsBlockTrans / 2.0 ) );
    G4double const dPhi = atan2( ta, r );

    // The r also adjusts because the detector blocks are flat (I assume)
    G4double const dR = ta * sin( dPhi ) / 2.0;

    // Set the translation
    G4ThreeVector position;
    position.setRhoPhiZ( r + dR, phi + dPhi, z + dZ );
    m_positions.push_back( position );

    // Set the rotation
    G4RotationMatrix * rotation = new G4RotationMatrix();
    rotation->rotateZ( -phi );
    m_rotations.push_back( rotation );

    // Visual properties
    G4VisAttributes* vis = new G4VisAttributes();
    //vis->SetColor( 0.0, G4double( copyNo ) / G4double( nCopies ), 0.0, 1.0 );
    vis->SetColor( 0.0, G4double( rand() % nCopies ) / G4double( nCopies ), 0.0, 1.0 );
    m_visions.push_back( vis );
  }
}

void ExplorerParameterisationCrystals::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  if ( copyNo >= ( G4int )m_positions.size() )
  {
    G4cerr << "Unknown copyNo for ExplorerParameterisationCrystals: " << copyNo << G4endl;
    return;
  }

  // Return precalculated result
  physVol->SetTranslation( m_positions.at( copyNo ) );
  physVol->SetRotation( m_rotations.at( copyNo ) );
  physVol->GetLogicalVolume()->SetVisAttributes( m_visions.at( copyNo ) );
}
