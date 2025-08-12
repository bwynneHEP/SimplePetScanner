#include "SiemensQuadraDetector.h"
#include "SiemensQuadraParameterisationCrystals.h"
#include "SiemensQuadraParameterisationBlocks.h"
#include "SiemensQuadraParameterisationPanels.h"
#include "CrystalMaterial.h"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

G4VPhysicalVolume* SiemensQuadraDetector::Construct( std::string Name, G4LogicalVolume* worldLV, std::string Mode, EnergyCounter * Counter, G4double LengthMM, std::string Material )
{
  // Default length
  if ( LengthMM <= 0.0 )
  {
    LengthMM = 1024.0;
  }

  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* air = nistManager->FindOrBuildMaterial( "G4_AIR" );
  G4Material* crystal = CrystalMaterial::GetCrystalMaterial(Material, "Siemens");
 
  std::cout << "Selected detector material: " << crystal->GetName() << std::endl;

  G4int nRings = NRingsInLength( LengthMM );
  if ( nRings == 32 )
  {
    std::cout << "Siemens Quadra detector with nRings: " << nRings << std::endl;
  }
  else
  {
    std::cout << "Siemens Quadra detector variant with nRings: " << nRings << std::endl;
  }

  // Panels of blocks, contiguous in axial direction
  G4double const panelAxial = LengthForNRings( nRings ) / 2.0; // half-length as always

  // Cylindrical envelope to contain whole detector
  // (non-physical, allows use of parameterised detector crystals)
  G4double const envelopeInnerRadius = 38.0 * cm;
  G4double const envelopeOuterRadius = 50.0 * cm;
  G4double const envelopeAxial = blockAxial * ( nRings + 1 );

  G4double y = crystalWidth;
  G4double z = crystalWidth;
  if ( Mode == "Block" || Mode == "BlockEnergy" )
  {
    y = blockTrans;
    z = blockAxial;
  }
  else if ( Mode == "Panel" )
  {
    y = blockTrans;
    z = panelAxial;
  }

  // ENVELOPE: Solid (cylinder)
  G4Tubs* envelopeS = new G4Tubs(
                 "Envelope",      // its name
                 envelopeInnerRadius, // inner radius, so it's a hollow tube
                 envelopeOuterRadius,   // outer radius
                 envelopeAxial,   // how long (z axis, remember it's a half-length)
                 0.0*deg,         // starting angle
                 360.0*deg );     // ending angle (i.e. it's a full circle)

  // ENVELOPE: Logical volume (how to treat it)
  G4LogicalVolume* envelopeLV = new G4LogicalVolume(
                 envelopeS,       // its solid
                 air,             // its material
                 "Envelope",      // its name
                 0, 0, 0 );       // Modifiers we don't use

  // ENVELOPE: Physical volume (where is it)
  /*G4VPhysicalVolume* envelopePV =*/ new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(0.0, 0.0, 0.0), // in the centre
                 envelopeLV,      // its logical volume
                 "Envelope",      // its name
                 worldLV,         // its mother volume
                 false,           // no boolean operations
                 0,               // copy number
                 true );          // checking overlaps

  // DETECTOR: the solid shape
  G4Box* detectorS = new G4Box(
                 Name,
                 crystalLength,
                 y,
                 z );

  // DETECTOR: Logical volume (how to treat it)
  G4LogicalVolume* detectorLV = new G4LogicalVolume(
                 detectorS,         // its solid
                 crystal,           // its material
                 Name,              // its name
                 0, 0, 0 );         // Modifiers we don't use

  // DETECTOR: Physical volume, parameterised to copy, rotate and translate the crystals
  G4int blocksPerRing = 38;
  if ( Mode == "Crystal" )
  {
    G4VPVParameterisation* detectorParam = new SiemensQuadraParameterisationCrystals( 200*blocksPerRing*nRings, Counter );
    return new G4PVParameterised( Name, detectorLV, envelopeLV, kUndefined, 200*blocksPerRing*nRings, detectorParam );
  }
  else if ( Mode == "CrystalNoIDs" )
  {
    G4VPVParameterisation* detectorParam = new SiemensQuadraParameterisationCrystals( 200*blocksPerRing*nRings, nullptr );
    return new G4PVParameterised( Name, detectorLV, envelopeLV, kUndefined, 200*blocksPerRing*nRings, detectorParam );
  }
  else if ( Mode == "Block" )
  {
    G4VPVParameterisation* detectorParam = new SiemensQuadraParameterisationBlocks( blocksPerRing*nRings, nullptr );
    return new G4PVParameterised( Name, detectorLV, envelopeLV, kUndefined, blocksPerRing*nRings, detectorParam );
  }
  else if ( Mode == "BlockEnergy" )
  {
    G4VPVParameterisation* detectorParam = new SiemensQuadraParameterisationBlocks( blocksPerRing*nRings, Counter );
    return new G4PVParameterised( Name, detectorLV, envelopeLV, kUndefined, blocksPerRing*nRings, detectorParam );
  }
  else if ( Mode == "Panel" )
  {
    G4VPVParameterisation* detectorParam = new SiemensQuadraParameterisationPanels( blocksPerRing );
    return new G4PVParameterised( Name, detectorLV, envelopeLV, kUndefined, blocksPerRing, detectorParam );
  }
  else
  {
    std::cerr << "Unrecognised Siemens Quadra detector mode: " << Mode << std::endl;
    exit(1);
  }
}

G4int SiemensQuadraDetector::NRingsInLength( G4double const LengthMM )
{
  return ceil( LengthMM * mm / ( blockAxial * 2.0 ) );
}
G4double SiemensQuadraDetector::LengthForNRings( G4int const NRings )
{
  return blockAxial * 2.0 * NRings;
}
