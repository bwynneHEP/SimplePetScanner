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
#include "G4PVParameterised.hh"

G4VPhysicalVolume* ExplorerDetector::Construct( std::string Name, G4LogicalVolume* worldLV, std::string Mode, G4double LengthMM, std::string Material )
{
  // Default length
  if ( LengthMM <= 0.0 )
  {
    LengthMM = 1872.92;//1940.0;
  }

  // Materials
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* air = nistManager->FindOrBuildMaterial( "G4_AIR" );
  G4bool isotopes = false;

  G4Material* crystal = nullptr;

  G4Element* O  = nistManager->FindOrBuildElement( "O" , isotopes );
  G4Element* Si = nistManager->FindOrBuildElement( "Si", isotopes );
  G4Element* Lu = nistManager->FindOrBuildElement( "Lu", isotopes );
  G4Material* LSO = new G4Material( "Lu2SiO5", 7.4*g/cm3, 3 );
  LSO->AddElement( Lu, 2 );
  LSO->AddElement( Si, 1 );
  LSO->AddElement( O , 5 );

  G4Element* Y  = nistManager->FindOrBuildElement( "Y" , isotopes );
  G4Material* LYSO = new G4Material( "LYSO", 7.1*g/cm3, 4 );
  LYSO->AddElement( Lu, 71.447 * perCent );
  LYSO->AddElement( Y,  4.034  * perCent );
  LYSO->AddElement( Si, 6.371  * perCent );
  LYSO->AddElement( O,  18.148 * perCent );

  //Hanna: material available via nistManager, no need to define
  G4Material* NaI = nistManager->FindOrBuildMaterial("G4_SODIUM_IODIDE");

  //Hanna: material available via nistManager, no need to define
  G4Material* BGO = nistManager->FindOrBuildMaterial("G4_BGO");

  //Hanna: material available via nistManager, no need to define
  G4Material* CsF = nistManager->FindOrBuildMaterial("G4_CESIUM_FLUORIDE");

  std::map<std::string, G4Material*> materialMap = {
    {"LSO", LSO},
    {"LYSO", LYSO},
    {"NaI", NaI},
    {"BGO", BGO}, 
    {"CsF", CsF}
  };

  if (auto search = materialMap.find(Material); search != materialMap.end())
        crystal = search->second;
  else if (Material.empty())
    crystal = LYSO;
  else {
    std::cerr << "Unrecognised detector material: " << Material << std::endl;
    exit(1);
  }

  std::cout << "Selected detector material: " << crystal->GetName() << std::endl;

  G4int nRings = NRingsInLength( LengthMM );
  if ( nRings == 8 )
  {
    std::cout << "Explorer detector with nRings: " << nRings << std::endl;
  }
  else
  {
    std::cout << "Explorer detector variant with nRings: " << nRings << std::endl;
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
                 crystal,           // its material
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

G4int ExplorerDetector::NRingsInLength( G4double const LengthMM )
{
  return ceil( ( ( LengthMM * mm ) + blockOffset ) / ( ( blockAxial * 2.0 ) + blockOffset ) );
}
G4double ExplorerDetector::LengthForNRings( G4int const NRings )
{
  return ( blockAxial * 2.0 * NRings ) + ( blockOffset * ( NRings - 1 ) );
}
