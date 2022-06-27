#include "SiemensQuadraDetector.h"
#include "SiemensQuadraParameterisation.h"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVParameterised.hh"

#include "G4RotationMatrix.hh"

G4VPhysicalVolume* SiemensQuadraDetector::Construct( std::string Name, G4LogicalVolume* worldLV )
{
  // Materials
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool isotopes = false;

  // LSO
  G4Element* O  = nistManager->FindOrBuildElement( "O" , isotopes );
  G4Element* Si = nistManager->FindOrBuildElement( "Si", isotopes );
  G4Element* Lu = nistManager->FindOrBuildElement( "Lu", isotopes );

  G4Material* LSO = new G4Material( "Lu2SiO5", 7.4*g/cm3, 3 );
  LSO->AddElement( Lu, 2 );
  LSO->AddElement( Si, 1 );
  LSO->AddElement( O , 5 );

  // Definitions of Solids, Logical Volumes, Physical Volumes
  G4double crystalWidth = 3.2*mm / 2.0; // half because it's measured from middle to face
  G4double crystalLength = 20.0*mm / 2.0;
  G4double blockAxial = crystalWidth * 10.0;
  G4double blockTrans = crystalWidth * 20.0;

  // DETECTOR: Single crystal (square prism)
  G4Box* detectorS = new G4Box(
                 Name,
                 crystalLength,
                 crystalWidth,
                 crystalWidth );

  // DETECTOR: 5x5 mini-blocks of crystals, 2x4 blocks
  /*G4Box* detectorS = new G4Box(
                 Name,
                 crystalLength,
                 blockTrans,
                 blockAxial );*/

  // DETECTOR: Logical volume (how to treat it)
  G4LogicalVolume* detectorLV = new G4LogicalVolume(
                 detectorS,         // its solid
                 LSO,               // its material
                 Name,              // its name
                 0, 0, 0 );         // Modifiers we don't use

  // DETECTOR: Physical volume, parameterised to copy, rotate and translate the crystals
  G4VPVParameterisation* detectorParam = new SiemensQuadraParameterisation();
  G4VPhysicalVolume* detectorPV = new G4PVParameterised( Name, detectorLV, worldLV, kUndefined, 243200, detectorParam ); //crystals
  //G4VPhysicalVolume* detectorPV = new G4PVParameterised( Name, detectorLV, worldLV, kZAxis, 1216, detectorParam ); //blocks

  return detectorPV;
}
