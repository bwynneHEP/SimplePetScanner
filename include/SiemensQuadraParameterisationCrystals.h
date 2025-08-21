// This class describes how the detector crystals are laid out in space

#ifndef SiemensQuadraParameterisationCrystals_h
#define SiemensQuadraParameterisationCrystals_h 1

#include "EnergyCounter.h"
#include "DetectorGeometryWriter.h"

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

class SiemensQuadraParameterisationCrystals : public G4VPVParameterisation
{
  public:
    SiemensQuadraParameterisationCrystals( G4int nCopies, EnergyCounter * Counter, DetectorGeometryData * DetectorData );
    ~SiemensQuadraParameterisationCrystals(){};

    void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const override;

  private:
    std::vector< G4ThreeVector > m_positions;
    std::vector< G4RotationMatrix* > m_rotations;
    std::vector< G4VisAttributes* > m_visions;
    std::vector< std::vector< int > > m_geometryIDs;
    EnergyCounter * m_counter = nullptr;
};

#endif
