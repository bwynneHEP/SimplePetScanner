// Siemens Quadra detector (approximate)
// Properties as listed in Table 1: https://jnm.snmjournals.org/content/early/2021/07/22/jnumed.121.261972

#ifndef SiemensQuadraDetector_h
#define SiemensQuadraDetector_h 1

#include "G4VPhysicalVolume.hh"

class SiemensQuadraDetector
{
  public:
    static G4VPhysicalVolume* Construct( std::string Name, G4LogicalVolume* worldLV );
};

#endif
