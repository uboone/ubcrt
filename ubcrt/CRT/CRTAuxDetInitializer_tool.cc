#include "larcorealg/Geometry/AuxDetReadoutGeom.h"

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include <cstdint>
#include <sstream>

using std::uint32_t;

namespace {

  class CRTAuxDetInitializer : public geo::AuxDetInitializer {
  public:
    explicit CRTAuxDetInitializer(fhicl::ParameterSet const& pset)
      : fNModules{pset.get<uint32_t>("NModules", 76)}
      , fNStripsPerModule{pset.get<uint32_t>("NStripsPerModule", 16)}
    {}

  private:
    geo::AuxDetReadoutInitializers initialize(std::vector<geo::AuxDetGeo> const& ads) const override
    {
      /**
       * Sets up the following variables:
       *
       * fADGeoToName: maps the channel number to name
       *
       * fNameToADGeo: maps the name to channel number. This seems degenerate.
       *
       * fADGeoToChannelAndSV: largely unused since the geometry has 1 sensitive volume per detector
       *
       * \todo: The two for loops can be exchanged for a single iterator pattern.
       *        This will make the code more extensible.
       **/

      geo::AuxDetReadoutInitializers result;
      uint32_t index=0;
      for(uint32_t mod=0; mod < fNModules; ++mod) {
        for(uint32_t strip=0; strip < fNStripsPerModule; ++strip) {
          std::ostringstream stream;
          stream << "volAuxDet_Module_" << mod << "_strip_" << strip;
          result.ADGeoToName[index] = stream.str();
          result.NameToADGeo[stream.str()] = index;
          result.ADGeoToChannelAndSV[index].emplace_back(index, 0);
          ++index;
        }
      }
      return result;
    }

    uint32_t fNModules;
    uint32_t fNStripsPerModule;
  };

}

DEFINE_ART_CLASS_TOOL(CRTAuxDetInitializer)
