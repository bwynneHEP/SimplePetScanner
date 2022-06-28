#include "SiemensQuadraDetector.h"
#include "SiemensQuadraParameterisationCrystals.h"
#include "SiemensQuadraParameterisationBlocks.h"
#include "SiemensQuadraParameterisationPanels.h"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVParameterised.hh"

#include "G4RotationMatrix.hh"

G4VPhysicalVolume* SiemensQuadraDetector::Construct( std::string Name, G4LogicalVolume* worldLV, std::string Mode )
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

  // Single crystal (square prism)
  G4double const crystalWidth = 3.2*mm / 2.0; // half because it's measured from middle to face
  G4double const crystalLength = 20.0*mm / 2.0;

  // 5x5 mini-blocks of crystals, 2x4 blocks
  G4double const blockAxial = crystalWidth * 10.0;
  G4double const blockTrans = crystalWidth * 20.0;

  // Panels of blocks, contiguous in axial direction
  G4double const panelAxial = blockAxial * 32.0;

  G4double y = crystalWidth;
  G4double z = crystalWidth;
  if ( Mode == "block" )
  {
    y = blockTrans;
    z = blockAxial;
  }
  else if ( Mode == "panel" )
  {
    y = blockTrans;
    z = panelAxial;
  }

  // DETECTOR: the solid shape
  G4Box* detectorS = new G4Box(
                 Name,
                 crystalLength,
                 y,
                 z );

  // DETECTOR: Logical volume (how to treat it)
  G4LogicalVolume* detectorLV = new G4LogicalVolume(
                 detectorS,         // its solid
                 LSO,               // its material
                 Name,              // its name
                 0, 0, 0 );         // Modifiers we don't use

  // DETECTOR: Physical volume, parameterised to copy, rotate and translate the crystals
  if ( Mode == "crystal" )
  {
    G4VPVParameterisation* detectorParam = new SiemensQuadraParameterisationCrystals();
    return new G4PVParameterised( Name, detectorLV, worldLV, kUndefined, 243200, detectorParam );
  }
  else if ( Mode == "block" )
  {
    G4VPVParameterisation* detectorParam = new SiemensQuadraParameterisationBlocks();
    return new G4PVParameterised( Name, detectorLV, worldLV, kUndefined, 1216, detectorParam );
  }
  else if ( Mode == "panel" )
  {
    G4VPVParameterisation* detectorParam = new SiemensQuadraParameterisationPanels();
    return new G4PVParameterised( Name, detectorLV, worldLV, kUndefined, 38, detectorParam );
  }
  else
  {
    G4cerr << "Unrecognised Siemens Quadra detector mode: " << Mode << G4endl;
    return nullptr;
  }
}
