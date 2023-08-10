#include "CrystalMaterial.h"

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

namespace CrystalMaterial{

G4Material* GetCrystalMaterial(const std::string& Material){
    // Materials
  G4NistManager* nistManager = G4NistManager::Instance();
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

  G4Material* CsI = nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");

  G4Material* BaF2 = nistManager->FindOrBuildMaterial("G4_BARIUM_FLUORIDE");

  std::map<std::string, G4Material*> materialMap = {
    {"LSO", LSO},
    {"LYSO", LYSO},
    {"NaI", NaI},
    {"BGO", BGO}, 
    {"CsF", CsF},
    {"CsI", CsI},
    {"BaF2", BaF2}
  };

  if (auto search = materialMap.find(Material); search != materialMap.end())
        crystal = search->second;
  else if (Material == "")
    crystal = LSO;
  else {
    std::cerr << "Unrecognised detector material: " << Material << std::endl;
    exit(1);
  }
  return crystal;
}
}