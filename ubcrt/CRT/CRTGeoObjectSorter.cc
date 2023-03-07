#include "ubcrt/CRT/CRTGeoObjectSorter.hh"
#include "larcorealg/Geometry/AuxDetGeo.h"

#include "fhiclcpp/ParameterSet.h"

#include <sstream>
#include <string>

namespace crt {

  CRTGeoObjectSorter::CRTGeoObjectSorter(
      fhicl::ParameterSet const& pSet):
      fNModules(pSet.get<uint32_t>("NModules", 76)),
      fNStripsPerModule(pSet.get<uint32_t>("NStripsPerModule", 16))
  {}
  
  bool CRTGeoObjectSorter::compareAuxDets(geo::AuxDetGeo const& ad1, geo::AuxDetGeo const& ad2) const
  {
    std::string name1 = ad1.Name();
    std::string name2 = ad2.Name();
    uint32_t mod1=0, mod2=0, strip1=0, strip2=0;
    for(uint32_t i=0; i<fNModules;++i){
      std::ostringstream stream;
      stream<<"Module_"<<i;
      if(name1.find(stream.str())!= std::string::npos) mod1 = i;
      if(name2.find(stream.str())!= std::string::npos) mod2 = i;
    }
    if(mod1>mod2) return true;
    if(mod2>mod1) return false;
    for(uint32_t i=0; i<fNStripsPerModule;++i){
      std::ostringstream stream;
      stream<<"strip_"<<i;
      if(name1.find(stream.str())!= std::string::npos) strip1 = i;
      if(name2.find(stream.str())!= std::string::npos) strip2 = i;
    }
    return strip1>strip2;
  }

  bool CRTGeoObjectSorter::compareAuxDetSensitives(geo::AuxDetSensitiveGeo const& ads1,
                                                   geo::AuxDetSensitiveGeo const& ads2) const
  {
    /// There is 1 SV per AD so the return value is irrelevant.
    return false;
  }

}  /// namespace crt
