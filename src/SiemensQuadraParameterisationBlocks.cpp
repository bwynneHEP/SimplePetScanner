#include "SiemensQuadraParameterisationBlocks.h"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

SiemensQuadraParameterisationBlocks::SiemensQuadraParameterisationBlocks( G4int nCopies, EnergyCounter * Counter ) :
  m_counter( Counter )
{
  // Precalculating everything avoids a memory leak
  m_positions.reserve( nCopies );
  m_rotations.reserve( nCopies );
  m_visions.reserve( nCopies );

  for ( G4int copyNo = 0; copyNo < nCopies; ++copyNo )
  {
    // 32 rings in the detector
    G4int const blocksPerRing = 38;
    G4int const ring = floor( copyNo / blocksPerRing );
    G4int const nRings = ceil( nCopies / blocksPerRing );
    G4int const inRing = copyNo % blocksPerRing;

    // mini-blocks are 5x5 crystals, arranged into 2x4 blocks (2 in the axial direction I think)
    // blocks are therefore 10x20
    G4int const crystalsBlockAxial = 10;

    // Phi position is block within ring
    G4double const deltaPhi = 360.0 * deg / G4double( blocksPerRing );
    G4double phi = deltaPhi * G4double( inRing );

    // Z position is ring itself
    G4double const crystalWidth = 3.2 * mm;
    G4double const ringWidth = crystalWidth * crystalsBlockAxial;
    //G4double const z = ( G4double( ring ) - G4double( nRings - 1 ) / 2.0 ) * ( ringWidth + 1*mm ); // +1mm for easy view
    G4double const z = ( G4double( ring ) - G4double( nRings - 1 ) / 2.0 ) * ringWidth; // offset to get it centred


    G4double const r = 41.0 * cm; // 82cm "Detector ring diameter"

    // Set the translation
    G4ThreeVector position;
    position.setRhoPhiZ( r, phi, z );
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

void SiemensQuadraParameterisationBlocks::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  if ( copyNo >= ( G4int )m_positions.size() )
  {
    G4cerr << "Unknown copyNo for SiemensQuadraParameterisationBlocks: " << copyNo << G4endl;
    return;
  }

  // Just for fun, make the crystal colour scale with total energy
  if ( m_counter )
  {
    m_visions.at( copyNo )->SetColor( 0.0, m_counter->GetEFraction( copyNo ), 0.0, 1.0 );
  }

  // Return precalculated result
  physVol->SetTranslation( m_positions.at( copyNo ) );
  physVol->SetRotation( m_rotations.at( copyNo ) );
  physVol->GetLogicalVolume()->SetVisAttributes( m_visions.at( copyNo ) );
}
