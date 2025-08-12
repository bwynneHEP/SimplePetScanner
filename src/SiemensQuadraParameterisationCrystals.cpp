#include "SiemensQuadraParameterisationCrystals.h"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

SiemensQuadraParameterisationCrystals::SiemensQuadraParameterisationCrystals( G4int nCopies, EnergyCounter * Counter )
{
  // Precalculating everything avoids a memory leak
  m_positions.reserve( nCopies );
  m_rotations.reserve( nCopies );
  m_visions.reserve( nCopies );
  m_geometryIDs.reserve( nCopies );
  m_counter = Counter;

  for ( G4int copyNo = 0; copyNo < nCopies; ++copyNo )
  {
    // 32 rings in the detector
    G4int const crystalsPerRing = 7600;
    G4int const blocksPerRing = 38;
    G4int const ring = floor( copyNo / crystalsPerRing );
    G4int const nRings = ceil( nCopies / crystalsPerRing );
    G4int const inRing = copyNo % crystalsPerRing;

    // 38 detector blocks per ring
    G4int const crystalsPerBlock = 200;
    G4int const block = floor( inRing / crystalsPerBlock );
    G4int const inBlock = inRing % crystalsPerBlock;

    // For GATE compatibility, I think rsector is like block, and ring is like module or submodule. Crystals in block is crystalID
    m_geometryIDs.emplace_back( std::initializer_list<int>{ring, block, inBlock} );

    // mini-blocks are 5x5 crystals, arranged into 2x4 blocks (2 in the axial direction I think)
    // blocks are therefore 10x20
    G4int const crystalsBlockAxial = 10;
    G4int const crystalsBlockTrans = 20;
    //G4int const blockTrans = floor( inBlock / crystalsBlockAxial );
    //G4int const blockAxial = inBlock % crystalsBlockAxial;
    G4int const blockAxial = floor( inBlock / crystalsBlockTrans ); // ideally this swap will match GATE
    G4int const blockTrans = inBlock % crystalsBlockTrans;

    // Phi position is block within ring
    G4double const deltaPhi = 360.0 * deg / G4double( blocksPerRing );
    G4double const phi = deltaPhi * G4double( block );

    // Z position is ring itself
    G4double const crystalWidth = 3.2 * mm;
    G4double const ringWidth = crystalWidth * crystalsBlockAxial;
    //G4double const z = ( G4double( ring ) - G4double( nRings - 1 ) / 2.0 ) * ( ringWidth + 1*mm ); // +1mm for easy view
    G4double const z = ( G4double( ring ) - G4double( nRings - 1 ) / 2.0 ) * ringWidth; // offset to get it centred

    // Adjust Z position for the block axial
    //G4double const dZ = ( crystalWidth + 1*mm ) * ( blockAxial - ( crystalsBlockAxial / 2.0 ) ); // +1mm for easy view
    G4double const dZ = crystalWidth * ( blockAxial - ( crystalsBlockAxial / 2.0 ) );

    // Adjust phi position for block transaxial
    G4double const r = 41.0 * cm; // 82cm "Detector ring diameter"
    //G4double const ta = ( crystalWidth + 1*mm ) * ( blockTrans - ( crystalsBlockTrans / 2.0 ) ); // +1mm for easy view
    G4double const ta = crystalWidth * ( blockTrans - ( crystalsBlockTrans / 2.0 ) );
    G4double const dPhi = atan2( ta, r );

    // The r also adjusts because the detector blocks are flat (I assume)
    G4double const dR = ta * sin( dPhi ) / 2.0;

    // Set the translation
    G4ThreeVector position;
    //position.setRhoPhiZ( r + dR, phi + dPhi, z + dZ );
    position.setRhoPhiZ( r + dR, phi + dPhi - (90*deg), z + dZ ); // added rotation to match GATE
    m_positions.push_back( position );

    // Set the rotation
    G4RotationMatrix * rotation = new G4RotationMatrix();
    rotation->rotateZ( -position.getPhi() );
    m_rotations.push_back( rotation );

    // Visual properties
    G4VisAttributes* vis = new G4VisAttributes();
    //vis->SetColor( 0.0, G4double( copyNo ) / G4double( nCopies ), 0.0, 1.0 ); // uniform
    vis->SetColor( 0.0, G4double( rand() % nCopies ) / G4double( nCopies ), 0.0, 1.0 ); // random speckle
    //vis->SetColor( (float)inBlock / (float)crystalsPerBlock, (float)block / (float)blocksPerRing, (float)ring / (float)nRings, 1.0 ); // GATE debug
    m_visions.push_back( vis );
  }

  // Push geometry ID information to the SD for output
  if ( m_counter ) m_counter->SetGeometryIDs( m_geometryIDs );
}

void SiemensQuadraParameterisationCrystals::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  if ( copyNo >= ( G4int )m_positions.size() )
  {
    G4cerr << "Unknown copyNo for SiemensQuadraParameterisationCrystals: " << copyNo << G4endl;
    return;
  }

  // Return precalculated result
  physVol->SetTranslation( m_positions.at( copyNo ) );
  physVol->SetRotation( m_rotations.at( copyNo ) );
  physVol->GetLogicalVolume()->SetVisAttributes( m_visions.at( copyNo ) );
}
