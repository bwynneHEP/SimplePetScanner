#ifndef DetectorGeometryWriter_h
#define DetectorGeometryWriter_h 1

#include <string>

struct DetectorGeometryData
{
  int nRings = 0;
};

class DetectorGeometryWriter
{
  public:
    static void WriteSTIRheader( std::string const& fileName, DetectorGeometryData const& inputData );
};

#endif
