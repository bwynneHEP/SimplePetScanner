#ifndef DetectorGeometryWriter_h
#define DetectorGeometryWriter_h 1

#include <string>

struct DetectorGeometryData
{
  int nRings = 0;
  float ringInnerDiameter = 0.0;
  float ringGap = 0.0;

  int blocksPerRing = 0;
  
  int crystalsAxial = 0;
  int crystalsTrans = 0;

  float crystalRadialSize = 0.0;
  float crystalAxialSize = 0.0;
  float crystalTransSize = 0.0;
};

class DetectorGeometryWriter
{
  public:
    static void WriteSTIRheader( std::string const& fileName, DetectorGeometryData const& inputData );
};

#endif
