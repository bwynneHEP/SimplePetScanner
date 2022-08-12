#include "ExplorerDetector.h"
#include "ExplorerParameterisationCrystals.h"
#include "ExplorerParameterisationBlocks.h"
#include "ExplorerParameterisationPanels.h"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVParameterised.hh"

#include "G4RotationMatrix.hh"

G4VPhysicalVolume* ExplorerDetector::Construct( std::string Name, G4LogicalVolume* worldLV, std::string Mode, G4double LengthMM, std::string )
{
  // Default length
  if ( LengthMM <= 0.0 )
  {
    LengthMM = 1850.0;//1940.0;
  }

  // Materials
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* air = nistManager->FindOrBuildMaterial( "G4_AIR" );
  G4bool isotopes = false;

  // LYSO
  // exact composition for EXPLORER from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6354919/
  G4Element* O  = nistManager->FindOrBuildElement( "O" , isotopes );
  G4Element* Si = nistManager->FindOrBuildElement( "Si", isotopes );
  G4Element* Lu = nistManager->FindOrBuildElement( "Lu", isotopes );
  G4Element* Y  = nistManager->FindOrBuildElement( "Y" , isotopes );

  G4Material* LYSO = new G4Material( "LYSO", 7.1*g/cm3, 5 );
  LYSO->AddElement( Lu, 71.447 * perCent );
  LYSO->AddElement( Y,  4.034  * perCent );
  LYSO->AddElement( Si, 6.371  * perCent );
  LYSO->AddElement( O,  18.148 * perCent );

  // Single crystal (square prism)
  G4double const crystalWidth = 2.76*mm / 2.0; // half because it's measured from middle to face
  G4double const crystalLength = 18.1*mm / 2.0;

  // 6x7 mini-blocks of crystals, 14x5 blocks
  G4double const blockAxial = crystalWidth * 84.0;
  G4double const blockTrans = crystalWidth * 35.0;

  // Work out how many rings
  G4int nRings = ceil( LengthMM*mm / ( blockAxial * 2.0 ) );
  if ( nRings == 8 )
  {
    std::cout << "Explorer detector with nRings: " << nRings << std::endl;
  }
  else
  {
    std::cout << "Explorer detector variant with nRings: " << nRings << std::endl;
  }

  // Panels of blocks, contiguous in axial direction
  G4double const blockOffset = 2.6*mm;
  G4double const panelAxial = ( blockAxial * nRings ) + ( blockOffset * ( nRings - 1 ) );

  // Cylindrical envelope to contain whole detector
  // (non-physical, allows use of parameterised detector crystals)
  G4double const envelopeInnerRadius = 38.0 * cm;
  G4double const envelopeOuterRadius = 50.0 * cm;
  G4double const envelopeAxial = blockAxial * ( nRings + 1 );

  G4double y = crystalWidth;
  G4double z = crystalWidth;
  if ( Mode == "Block" )
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
                 LYSO,              // its material
                 Name,              // its name
                 0, 0, 0 );         // Modifiers we don't use

  // DETECTOR: Physical volume, parameterised to copy, rotate and translate the crystals
  G4int blocksPerRing = 24;
  if ( Mode == "Crystal" )
  {
    G4VPVParameterisation* detectorParam = new ExplorerParameterisationCrystals( 2940*blocksPerRing*nRings );
    return new G4PVParameterised( Name, detectorLV, envelopeLV, kUndefined, 2940*blocksPerRing*nRings, detectorParam );
  }
  else if ( Mode == "Block" )
  {
    G4VPVParameterisation* detectorParam = new ExplorerParameterisationBlocks( blocksPerRing*nRings );
    return new G4PVParameterised( Name, detectorLV, envelopeLV, kUndefined, blocksPerRing*nRings, detectorParam );
  }
  else if ( Mode == "Panel" )
  {
    G4VPVParameterisation* detectorParam = new ExplorerParameterisationPanels( blocksPerRing );
    return new G4PVParameterised( Name, detectorLV, envelopeLV, kUndefined, blocksPerRing, detectorParam );
  }
  else
  {
    std::cerr << "Unrecognised Explorer detector mode: " << Mode << std::endl;
    exit(1);
  }
}
