/**
 * \class CRTGeoObjectSorter
 *
 * \ingroup crt
 *
 * \brief Sorts AuDetGeos by module number and strip number
 *
 * In order to map the CRT panels, they need to be arranged by
 * module and strip number. This is accomplished by sorting on this
 * pattern.
 *
 * \author $Author: Kevin Wierman<kevin.wierman@pnnl.gov> 
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2016/12/12 $
 *
 * Contact: kevin.wierman@pnnl.gov
 *
 * Created on: Tuesday, December 13, 2016
 *
 */

#ifndef CRTGeoObjectSorter_hh_
#define CRTGeoObjectSorter_hh_

#include "larcorealg/Geometry/AuxDetGeoObjectSorter.h"
#include "fhiclcpp/fwd.h"

#include <cstdint>
#include <vector>

namespace crt {

  class CRTGeoObjectSorter : public geo::AuxDetGeoObjectSorter {
    std::uint32_t fNModules;
    std::uint32_t fNStripsPerModule;

  public:
    CRTGeoObjectSorter(fhicl::ParameterSet const& p);

    bool compareAuxDets(geo::AuxDetGeo const& ad1, geo::AuxDetGeo const& ad2) const override;
    bool compareAuxDetSensitives(geo::AuxDetSensitiveGeo const& ads1,
                                 geo::AuxDetSensitiveGeo const& ads2) const override;
  };
}

#endif  // CRTGeoObjectSorter_hh_