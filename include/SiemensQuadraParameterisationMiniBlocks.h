// This class describes how the detector crystals are laid out in space

#ifndef SiemensQuadraParameterisationMiniBlocks_h
#define SiemensQuadraParameterisationMiniBlocks_h 1

#include "EnergyCounter.h"

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

class SiemensQuadraParameterisationMiniBlocks : public G4VPVParameterisation
{
  public:
    SiemensQuadraParameterisationMiniBlocks( G4int nCopies, EnergyCounter * Counter );
    ~SiemensQuadraParameterisationMiniBlocks(){};

    void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const override;

  private:
    std::vector< G4ThreeVector > m_positions;
    std::vector< G4RotationMatrix* > m_rotations;
    std::vector< G4VisAttributes* > m_visions;
    EnergyCounter * m_counter;
};

#endif
