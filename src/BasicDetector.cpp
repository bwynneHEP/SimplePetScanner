#include "BasicDetector.h"
#include "BasicParameterisation.h"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVParameterised.hh"

G4VPhysicalVolume* BasicDetector::Construct( std::string Name, G4LogicalVolume* worldLV )
{
  // Materials
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* air = nistManager->FindOrBuildMaterial( "G4_AIR" );
  G4bool isotopes = false;

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

  // Cylindrical envelope to contain whole detector
  // (non-physical, allows use of parameterised detector crystals)
  G4double const envelopeInnerRadius = 39.0 * cm;
  G4double const envelopeOuterRadius = 61.0 * cm;
  G4double const envelopeAxial = 150.0 * cm;

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

  // DETECTOR: Single crystal (square prism)
  G4Box* detectorS = new G4Box(
                 Name,
                 detectorLength,
                 detectorWidth,
                 detectorWidth );

  // DETECTOR: Logical volume (how to treat it)
  G4LogicalVolume* detectorLV = new G4LogicalVolume(
                 detectorS,         // its solid
                 LYSO,              // its material
                 Name,              // its name
                 0, 0, 0 );         // Modifiers we don't use

  // DETECTOR: Physical volume, parameterised to copy, rotate and translate the crystals
  G4VPVParameterisation* detectorParam = new BasicParameterisation( 20 );
  G4VPhysicalVolume* detectorPV = new G4PVParameterised( Name, detectorLV, envelopeLV, kUndefined, 100, detectorParam );

  return detectorPV;
}
