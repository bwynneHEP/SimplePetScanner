#ifndef CrystalMaterial_h
#define CrystalMaterial_h 1

#include "G4Material.hh"

#include <string>

namespace CrystalMaterial{

G4Material* GetCrystalMaterial(const std::string& Material);

}

#endif