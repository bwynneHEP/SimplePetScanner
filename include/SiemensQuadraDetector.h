// Siemens Quadra detector (approximate)
// Properties as listed in Table 1: https://jnm.snmjournals.org/content/early/2021/07/22/jnumed.121.261972

#ifndef SiemensQuadraDetector_h
#define SiemensQuadraDetector_h 1

#include "EnergyCounter.h"
#include "DetectorGeometryWriter.h"

#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"

class SiemensQuadraDetector
{
  public:
    static G4VPhysicalVolume* Construct( std::string Name, G4LogicalVolume* worldLV, std::string Mode, EnergyCounter * Counter, DetectorGeometryData * DetectorData,
                                         G4double LengthMM=1024, std::string Material="LSO" );
    static G4int NRingsInLength( G4double const Length );
    static G4double LengthForNRings( G4int const NRings );

  private:
    // Single crystal (square prism)
    static G4double constexpr crystalWidth = 3.2*mm / 2.0; // half because it's measured from middle to face
    static G4double constexpr crystalLength = 20.0*mm / 2.0;

    // 5x5 mini-blocks of crystals, 2x4 blocks
    static G4double constexpr blockAxial = crystalWidth * 10.0;
    static G4double constexpr blockTrans = crystalWidth * 20.0;
};

#endif
