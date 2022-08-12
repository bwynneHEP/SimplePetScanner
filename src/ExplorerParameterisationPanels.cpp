#include "ExplorerParameterisationPanels.h"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

ExplorerParameterisationPanels::ExplorerParameterisationPanels( G4int nCopies )
{
  // Precalculating everything avoids a memory leak
  m_positions.reserve( nCopies );
  m_rotations.reserve( nCopies );
  m_visions.reserve( nCopies );

  for ( G4int copyNo = 0; copyNo < nCopies; ++copyNo )
  {
    // 32 rings in the detector
    G4int const blocksPerRing = 24;

    // Phi position is block within ring
    G4double const deltaPhi = 360.0 * deg / G4double( blocksPerRing );
    G4double const phi = deltaPhi * G4double( copyNo );

    G4double const z = 0.0; // panels go the full length

    G4double const r = 40.205 * cm; // 78.6 cm in diameter (detector face-to-face) + 18.1/2 crystal depth

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

void ExplorerParameterisationPanels::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  if ( copyNo >= ( G4int )m_positions.size() )
  {
    G4cerr << "Unknown copyNo for ExplorerParameterisationPanels: " << copyNo << G4endl;
    return;
  }

  // Return precalculated result
  physVol->SetTranslation( m_positions.at( copyNo ) );
  physVol->SetRotation( m_rotations.at( copyNo ) );
  physVol->GetLogicalVolume()->SetVisAttributes( m_visions.at( copyNo ) );
}
