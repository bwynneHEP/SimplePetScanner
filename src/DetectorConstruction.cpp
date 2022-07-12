#include "DetectorConstruction.h"
#include "EnergyCounter.h"
#include "BasicDetector.h"
#include "SiemensQuadraDetector.h"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::m_magneticFieldMessenger = 0;

DetectorConstruction::DetectorConstruction( DecayTimeFinderAction * decayTimeFinder, std::string detector )
  : G4VUserDetectorConstruction()
  , m_decayTimeFinder( decayTimeFinder )
  , m_detector( detector )
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
  G4Material* water = nistManager->FindOrBuildMaterial( "G4_WATER" );
  //G4Material* brain = nistManager->FindOrBuildMaterial("G4_BRAIN_ICRP");

  // Definitions of Solids, Logical Volumes, Physical Volumes
  G4double const worldAxial = 1.5*m;
  G4double const worldTrans = 1.0*m;
  G4double const phantomRadius = 20.3*cm / 2.0;
  G4double const phantomAxial = 35.0*cm;

  // WORLD: Solid (cube)
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent( worldAxial );
  G4Box* worldS = new G4Box(
                 "World",         // its name
                 worldTrans,
                 worldTrans,
                 worldAxial );   // its size (in half-lengths)

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

  // PHANTOM: Solid (cylinder)
  G4Tubs* phantomS = new G4Tubs(
                 "Phantom",       // its name
                 0.0,             // inner radius 0, so it's a solid cylinder (not a hollow tube)
                 phantomRadius,   // outer radius
                 phantomAxial,    // how long (z axis, remember it's a half-length)
                 0.0*deg,         // starting angle
                 360.0*deg );     // ending angle (i.e. it's a full circle)

  // PHANTOM: Logical volume (how to treat it)
  G4LogicalVolume* phantomLV = new G4LogicalVolume(
                 phantomS,        // its solid
                 water,           // its material
                 "Phantom",       // its name
                 0, 0, 0 );       // Modifiers we don't use

  // PHANTOM: Physical volume (where is it)
  /*G4VPhysicalVolume* phantomPV =*/ new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(0.0, 0.0, 0.0), // in the centre
                 phantomLV,       // its logical volume
                 "Phantom",       // its name
                 worldLV,         // its mother volume
                 false,           // no boolean operations
                 0,               // copy number
                 true );          // checking overlaps

  // DETECTOR: Physical volume, parameterised to copy, rotate and translate the crystals
  if ( m_detector == "SiemensCrystal" ) SiemensQuadraDetector::Construct( "Detector", worldLV, "crystal" );
  else if ( m_detector == "SiemensBlock" ) SiemensQuadraDetector::Construct( "Detector", worldLV, "block" );
  else if ( m_detector == "SiemensPanel" ) SiemensQuadraDetector::Construct( "Detector", worldLV, "panel" );
  else if ( m_detector == "Basic" ) BasicDetector::Construct( "Detector", worldLV );
  else G4cerr << "Unrecognised detector name: " << m_detector << G4endl;

  // DETECTOR: Warn if there's an overlap, but it's very slow (N^2 I think)
  //if ( detectorPV->CheckOverlaps() ) std::cerr << "WARNING: your simulated objects overlap" << std::endl;

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

  auto detector = new EnergyCounter( "Detector", m_decayTimeFinder );
  G4SDManager::GetSDMpointer()->AddNewDetector( detector );
  this->SetSensitiveDetector( "Detector", detector );
}
