#include "CrystalMaterial.h"

#include "G4NistManager.hh"
#include "G4Isotope.hh"
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

  //non-radioactive CsI
  G4Material* CsI = nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");

  //radioactive CsI with a fraction of Cs-137
  G4Isotope* Cs133 = new G4Isotope("Cs133", 55, 133, 132.905*g/mole);
  G4Isotope* Cs137 = new G4Isotope("Cs137", 55, 137, 136.907*g/mole);
  G4Element* eCs = new G4Element("Enriched Caesium", "Cs", 2);
  eCs->AddIsotope(Cs133, 99.*perCent);
  eCs->AddIsotope(Cs137, 1.*perCent);
  G4Element* I = nistManager->FindOrBuildElement( "I", isotopes );
  G4Material* eCsI = new G4Material("eCsI", 4.51*g/cm3, 2);
  eCsI->AddElement(eCs, 1);
  eCsI->AddElement(I, 1);
  // std::cout << *(G4Isotope::GetIsotopeTable()) << std::endl;
  // std::cout << *(G4Element::GetElementTable()) << std::endl;
  // std::cout << *(G4Material::GetMaterialTable()) << std::endl;

  G4Material* BaF2 = nistManager->FindOrBuildMaterial("G4_BARIUM_FLUORIDE");

  G4Material* CaF2 = nistManager->FindOrBuildMaterial("G4_CALCIUM_FLUORIDE");

  G4Element* Cd = nistManager->FindOrBuildElement( "Cd", isotopes );
  G4Element* W = nistManager->FindOrBuildElement( "W", isotopes );
  G4Material* CdWO4 = new G4Material( "CdWO4", 7.9*g/cm3, 3 );
  CdWO4->AddElement( Cd, 1 );
  CdWO4->AddElement( W, 1 );
  CdWO4->AddElement( O, 4 );

  G4Element* Zn = nistManager->FindOrBuildElement( "Zn", isotopes );
  G4Element* Te = nistManager->FindOrBuildElement( "Te", isotopes );
  G4Material* CZT = new G4Material( "CdZnTe", 5.8*g/cm3, 3 );
  CZT->AddElement( Cd, 1 );
  CZT->AddElement( Zn, 1 );
  CZT->AddElement( Te, 1 );

  //Perovskites
  G4Element* C  = nistManager->FindOrBuildElement( "C" , isotopes );
  G4Element* H = nistManager->FindOrBuildElement( "H", isotopes );
  G4Element* N = nistManager->FindOrBuildElement( "N", isotopes );
  G4Element* Pb  = nistManager->FindOrBuildElement( "Pb" , isotopes );
  G4Element* Br  = nistManager->FindOrBuildElement( "Br" , isotopes );
  G4Material* MAPbBr3 = new G4Material( "CH3NH3PbBr3", 3.83*g/cm3, 5 );
  MAPbBr3->AddElement( C, 2.51 * perCent );
  MAPbBr3->AddElement( H,  1.27  * perCent );
  MAPbBr3->AddElement( N, 2.92  * perCent );
  MAPbBr3->AddElement( Pb,  43.26 * perCent );
  MAPbBr3->AddElement(Br,  50.04 * perCent );

  G4Element* Cs = nistManager->FindOrBuildElement( "Cs" , isotopes );
  G4Element* Ag = nistManager->FindOrBuildElement( "Ag", isotopes );
  G4Element* Bi = nistManager->FindOrBuildElement( "Bi", isotopes );
  G4Material* Cs2AgBiBr6 = new G4Material( "Cs2AgBiBr6", 4.65*g/cm3, 4 );
  Cs2AgBiBr6->AddElement( Cs, 25.03 * perCent );
  Cs2AgBiBr6->AddElement( Ag,  10.16  * perCent );
  Cs2AgBiBr6->AddElement( Bi, 19.68  * perCent );
  Cs2AgBiBr6->AddElement( Br,  45.13 * perCent );

  G4Material* MAPbI3 = new G4Material( "CH3NH3PbI3", 4.0*g/cm3, 5 );
  MAPbI3->AddElement( C, 1 );
  MAPbI3->AddElement( H, 6 );
  MAPbI3->AddElement( N, 1 );
  MAPbI3->AddElement( Pb, 1 );
  MAPbI3->AddElement( I, 3 );

  G4Element* Ca = nistManager->FindOrBuildElement( "Ca", isotopes );
  G4Element* Ti = nistManager->FindOrBuildElement( "Ti", isotopes );
  G4Material* CaTiO3 = new G4Material( "CaTiO3", 0.986*g/cm3, 3 );
  CaTiO3->AddElement( Ca, 1 );
  CaTiO3->AddElement( Ti, 1 );
  CaTiO3->AddElement( O, 3 );

  G4Material* CsPbBr3 = new G4Material( "CsPbBr3", 5.8*g/cm3, 3 );
  CsPbBr3->AddElement( Cs, 1 );
  CsPbBr3->AddElement( Pb, 1 );
  CsPbBr3->AddElement( Br, 3 );

  G4Material* FAPbI3 = new G4Material( "(CH(NH2)2PbI3", 4.0*g/cm3, 5 ); //(CH(NH2)2PbI3
  FAPbI3->AddElement( C, 2 );
  FAPbI3->AddElement( H, 6 );
  FAPbI3->AddElement( N, 2 );
  FAPbI3->AddElement( Pb, 1 );
  FAPbI3->AddElement( I, 3 );

  std::map<std::string, G4Material*> materialMap = {
    {"LSO", LSO},
    {"LYSO", LYSO},
    {"NaI", NaI},
    {"BGO", BGO}, 
    {"CsF", CsF},
    {"CsI", CsI},
    {"eCsI", eCsI},
    {"BaF2", BaF2},
    {"CaF2", CaF2},
    {"CdWO4", CdWO4},
    {"CZT", CZT},
    {"MAPbBr3", MAPbBr3},
    {"Cs2AgBiBr6", Cs2AgBiBr6},
    {"MAPbI3", MAPbI3},
    {"CaTiO3", CaTiO3},
    {"CsPbBr3", CsPbBr3},
    {"FAPbI3", FAPbI3}
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