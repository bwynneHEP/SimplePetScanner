#include "DetectorConstruction.h"
#include "BasicDetector.h"
#include "SiemensQuadraDetector.h"
#include "ExplorerDetector.h"

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

DetectorConstruction::DetectorConstruction( DecayTimeFinderAction * decayTimeFinder, std::string detector, G4double detectorLength, G4double phantomLength,
                                            std::string outputFileName, std::string decayOutputFileName, std::string material )
  : G4VUserDetectorConstruction()
  , m_decayTimeFinder( decayTimeFinder )
  , m_detector( detector )
  , m_outputFileName( outputFileName )
  , m_decayOutputFileName( decayOutputFileName )
  , m_material( material )
  , m_detectorLength( detectorLength )
  , m_phantomLength( phantomLength )
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
  G4Material* polyeth = nistManager->FindOrBuildMaterial( "G4_POLYETHYLENE" );

  // Definitions of Solids, Logical Volumes, Physical Volumes
  G4double const worldAxial = 1.5*m;
  G4double const worldTrans = 1.0*m;
  //G4double const phantomRadius = 20.3*cm / 2.0;
  G4double const phantomRadius = 8.3*cm / 2.0;
  G4double phantomAxial = 35.0*cm;
  if ( m_phantomLength >= 0.0 ) phantomAxial = m_phantomLength * mm / 2.0; // half-lengths

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

  if ( phantomAxial > 0.0 ) {
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
                   polyeth,         // its material
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

		// Reservoir 1
    // PHANTOM: Solid (cylinder)
    G4Tubs* phantomR1S = new G4Tubs(
                   "PhantomR1",       // its name
                   0.0,             // inner radius 0, so it's a solid cylinder (not a hollow tube)
                   phantomRadius-(20*mm),   // outer radius
                   (phantomAxial/2.0)-(20*mm),    // how long (z axis, remember it's a half-length)
                   0.0*deg,         // starting angle
                   360.0*deg );     // ending angle (i.e. it's a full circle)

    // PHANTOM: Logical volume (how to treat it)
    G4LogicalVolume* phantomR1LV = new G4LogicalVolume(
                   phantomR1S,        // its solid
                   air,         // its material
                   "PhantomR1",       // its name
                   0, 0, 0 );       // Modifiers we don't use

    // PHANTOM: Physical volume (where is it)
    /*G4VPhysicalVolume* phantomPV =*/ new G4PVPlacement(
                   0,               // no rotation
                   G4ThreeVector(0.0, 0.0, phantomAxial/2.0), // in the centre
                   phantomR1LV,       // its logical volume
                   "PhantomR1",       // its name
                   phantomLV,         // its mother volume
                   false,           // no boolean operations
                   0,               // copy number
                   true );          // checking overlaps

    // Reservoir 2
    // PHANTOM: Solid (cylinder)
    G4Tubs* phantomR2S = new G4Tubs(
                   "PhantomR2",       // its name
                   0.0,             // inner radius 0, so it's a solid cylinder (not a hollow tube)
                   phantomRadius-(20*mm),   // outer radius
                   (phantomAxial/2.0)-(20*mm),    // how long (z axis, remember it's a half-length)
                   0.0*deg,         // starting angle
                   360.0*deg );     // ending angle (i.e. it's a full circle)

    // PHANTOM: Logical volume (how to treat it)
    G4LogicalVolume* phantomR2LV = new G4LogicalVolume(
                   phantomR2S,        // its solid
                   air,         // its material
                   "PhantomR2",       // its name
                   0, 0, 0 );       // Modifiers we don't use

    // PHANTOM: Physical volume (where is it)
    /*G4VPhysicalVolume* phantomPV =*/ new G4PVPlacement(
                   0,               // no rotation
                   G4ThreeVector(0.0, 0.0, -phantomAxial/2.0), // in the centre
                   phantomR2LV,       // its logical volume
                   "PhantomR2",       // its name
                   phantomLV,         // its mother volume
                   false,           // no boolean operations
                   0,               // copy number
                   true );          // checking overlaps

    // Holes
    for ( int xIndex = 0; xIndex < 5; ++xIndex ) {
    for ( int yIndex = 0; yIndex < 5; ++yIndex ) {
    std::string holeName = "Hole" + std::to_string( xIndex ) + "." + std::to_string( yIndex );
    //G4double xPos = (-7.5*mm)+(5.0*mm*(G4double)xIndex);
    //G4double yPos = (-7.5*mm)+(5.0*mm*(G4double)yIndex);
    G4double xPos = (-10.0*mm)+(5.0*mm*(G4double)xIndex);
    G4double yPos = (-10.0*mm)+(5.0*mm*(G4double)yIndex);
    G4Tubs* holeS = new G4Tubs(
                   holeName,       // its name
                   0.0,             // inner radius 0, so it's a solid cylinder (not a hollow tube)
                   2*mm,   // outer radius
                   20*mm,    // how long (z axis, remember it's a half-length)
                   0.0*deg,         // starting angle
                   360.0*deg );     // ending angle (i.e. it's a full circle)

    // PHANTOM: Logical volume (how to treat it)
    G4LogicalVolume* holeLV = new G4LogicalVolume(
                   holeS,        // its solid
                   air,         // its material
                   holeName,       // its name
                   0, 0, 0 );       // Modifiers we don't use

    // PHANTOM: Physical volume (where is it)
    /*G4VPhysicalVolume* phantomPV =*/ new G4PVPlacement(
                   0,               // no rotation
                   G4ThreeVector(xPos, yPos, 0.0), // in the centre
                   holeLV,       // its logical volume
                   holeName,       // its name
                   phantomLV,         // its mother volume
                   false,           // no boolean operations
                   0,               // copy number
                   true );          // checking overlaps
                   }}
  }

  // DETECTOR: separate class
  m_energyCounter = new EnergyCounter( "Detector", m_decayTimeFinder, m_outputFileName );
  if ( m_detector.substr( 0, 7 ) == "Siemens" )
  {
    SiemensQuadraDetector::Construct( "Detector", worldLV, m_detector.substr( 7 ), m_energyCounter, m_detectorLength, m_material );
  }
  else if ( m_detector.substr( 0, 8 ) == "Explorer" )
  {
    ExplorerDetector::Construct( "Detector", worldLV, m_detector.substr( 8 ), m_detectorLength, m_material );
  }
  else if ( m_detector == "Basic" )
  {
    BasicDetector::Construct( "Detector", worldLV );
  }
  else
  {
    std::cerr << "Unrecognised detector name: " << m_detector << std::endl;
    exit(1);
  }

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

  // Sensitive detector
  G4SDManager::GetSDMpointer()->AddNewDetector( m_energyCounter );
  this->SetSensitiveDetector( "Detector", m_energyCounter );
}
