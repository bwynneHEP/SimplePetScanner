// Simple cylindrical detector

#ifndef BasicDetector_h
#define BasicDetector_h 1

#include "G4VPhysicalVolume.hh"

class BasicDetector
{
  public:
    static G4VPhysicalVolume* Construct( std::string Name, G4LogicalVolume* worldLV );
};

#endif
