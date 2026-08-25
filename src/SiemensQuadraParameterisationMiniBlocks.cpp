#include "SiemensQuadraParameterisationMiniBlocks.h"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

SiemensQuadraParameterisationMiniBlocks::SiemensQuadraParameterisationMiniBlocks( G4int nCopies, EnergyCounter * Counter ) :
  m_counter( Counter )
{
  // Precalculating everything avoids a memory leak
  m_positions.reserve( nCopies );
  m_rotations.reserve( nCopies );
  m_visions.reserve( nCopies );

  for ( G4int copyNo = 0; copyNo < nCopies; ++copyNo )
  {
    // Halfring axial width is half of the ring and therefore half of the axial block width, i.e. mini block width
    // Transaxially, ring has 38 blocks, one block has 4 miniblocks, so there are 152 mini blocks in a ring 
    G4int const miniBlocksPerRingTot = 304;
    G4int const miniBlocksPerRingAxially = 2;
    G4int const miniBlocksPerRingTrans = 152;
    G4int const halfring = floor( copyNo / miniBlocksPerRingTrans ); 
    G4int const nRings = ceil( nCopies / miniBlocksPerRingTot );
    G4int const nHalfrings = nRings*2;
    G4int const inHalfring = copyNo % miniBlocksPerRingTrans;

    G4double const deltaPhi = 360.0 * deg / G4double( miniBlocksPerRingTrans );
    G4double phi = deltaPhi * G4double( inHalfring );

    // // Z position is ring itself
    G4double const crystalWidth = 3.2 * mm;
    G4double const halfringWidth = 5.0*crystalWidth;
    G4double const z = ( G4double( halfring ) - G4double( nHalfrings - 1 ) / 2.0 ) * halfringWidth; // offset to get it centred;

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
    vis->SetColor( 0.0, G4double( rand() % nCopies ) / G4double( nCopies ), 0.0, 1.0 );
    m_visions.push_back( vis );
  }
}

void SiemensQuadraParameterisationMiniBlocks::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  if ( copyNo >= ( G4int )m_positions.size() )
  {
    G4cerr << "Unknown copyNo for SiemensQuadraParameterisationMiniBlocks: " << copyNo << G4endl;
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
