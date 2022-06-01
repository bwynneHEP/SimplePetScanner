#include "DetectorConstruction.h"
#include "EnergyCounter.h"
#include "Parameterisation.h"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4PVParameterised.hh"

// Set number of detector layers
G4int const nLayers = 5;

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::m_magneticFieldMessenger = 0;

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
}

DetectorConstruction::~DetectorConstruction()
{
}

// Here we define the actual experiment that we want to perform
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Materials
  // http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/Appendix/materialNames.html
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* air = nistManager->FindOrBuildMaterial( "G4_AIR" );

  G4bool isotopes = false;

  // LSO
  /*
  G4Element* O  = nistManager->FindOrBuildElement( "O" , isotopes );
  G4Element* Si = nistManager->FindOrBuildElement( "Si", isotopes );
  G4Element* Lu = nistManager->FindOrBuildElement( "Lu", isotopes );

  G4Material* LSO = new G4Material( "Lu2SiO5", 7.4*g/cm3, 3 );
  LSO->AddElement( Lu, 2 );
  LSO->AddElement( Si, 1 );
  LSO->AddElement( O , 5 );
  */

  // LYSO
  G4Element* O  = nistManager->FindOrBuildElement( "O" , isotopes );
  G4Element* Si = nistManager->FindOrBuildElement( "Si", isotopes );
  G4Element* Lu = nistManager->FindOrBuildElement( "Lu", isotopes );
  G4Element* Ce = nistManager->FindOrBuildElement( "Ce", isotopes );
  G4Element* Y  = nistManager->FindOrBuildElement( "Y" , isotopes );

  G4Material* LYSO = new G4Material( "LYSO", 7.1*g/cm3, 5 );
  LYSO->AddElement( Lu, 71.43 * perCent );
  LYSO->AddElement( Y,  4.03  * perCent );
  LYSO->AddElement( Si, 6.37  * perCent );
  LYSO->AddElement( O,  18.14 * perCent );
  LYSO->AddElement( Ce, 0.02  * perCent );

  // Definitions of Solids, Logical Volumes, Physical Volumes
  G4double detectorWidth = 5.0*cm;
  G4double detectorLength = 10.0*cm;
  G4double worldLength = 250.0*cm;

  // WORLD: Solid (cube)
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent( worldLength );
  G4Box* worldS = new G4Box(
                 "World",         // its name
                 worldLength,
                 worldLength,
                 worldLength );   // its size (in half-lengths)

  // WORLD: Logical volume (how to treat it)
  G4LogicalVolume* worldLV = new G4LogicalVolume(
                 worldS,          // its solid
                 air,             // its material
                 "World" );       // its name

  // WORLD: Physical volume (where is it)
  // Must place the World Physical volume unrotated at (0,0,0).
  G4VPhysicalVolume* worldPV = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(0.0, 0.0, 0.0), // in the centre
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother volume
                 false,           // no boolean operations
                 0,               // copy number
                 true );          // checking overlaps

  // DETECTOR: Single crystal (square prism)
  G4Box* detectorS = new G4Box(
                 "Detector",        // its name
                 detectorLength,
                 detectorWidth,
                 detectorWidth );

  // DETECTOR: Logical volume (how to treat it)
  G4LogicalVolume* detectorLV = new G4LogicalVolume(
                 detectorS,         // its solid
                 LYSO,              // its material
                 "Detector",        // its name
                 0, 0, 0 );         // Modifiers we don't use

  // DETECTOR: Physical volume, parameterised to copy, rotate and translate the crystals
  G4VPVParameterisation* detectorParam = new Parameterisation( 20 );
  G4VPhysicalVolume* detectorPV = new G4PVParameterised( "Detector", detectorLV, worldLV, kUndefined, 200, detectorParam );

  // DETECTOR: Warn if there's an overlap
  if ( detectorPV->CheckOverlaps() ) std::cerr << "WARNING: your simulated objects overlap" << std::endl;

  // Always return the physical world
  return worldPV;
}

// Set up the magnetic field
void DetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  m_magneticFieldMessenger = new G4GlobalMagFieldMessenger( fieldValue );

  // Register the field messenger for deleting
  G4AutoDelete::Register( m_magneticFieldMessenger );

  auto detector = new EnergyCounter( "Detector" );
  G4SDManager::GetSDMpointer()->AddNewDetector( detector );
  this->SetSensitiveDetector( "Detector", detector );
}
