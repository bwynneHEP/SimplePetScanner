// Explorer detector (approximate)
// Properties as listed in https://jnm.snmjournals.org/content/early/2020/10/23/jnumed.120.250597

#ifndef ExplorerDetector_h
#define ExplorerDetector_h 1

#include "G4VPhysicalVolume.hh"

class ExplorerDetector
{
  public:
    static G4VPhysicalVolume* Construct( std::string Name, G4LogicalVolume* worldLV, std::string Mode, G4double LengthMM=1940.0, std::string Material="LYSO" );
};

#endif
