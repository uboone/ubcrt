/**
 * \class CRTGeometryHelper
 *
 * \ingroup crt
 *
 * \brief Interface class for the crt channel map.
 *
 * See `geo::AuxDetExptGeoHelperInterface` for full explanation.
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

#ifndef CRTGeometryHelper_HH_
#define CRTGeometryHelper_HH_

// framework libraries
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/AuxDetExptGeoHelperInterface.h"
#include "larcorealg/Geometry/AuxDetChannelMapAlg.h"

#include <memory> //For std::shared_ptr


namespace crt
{

  class CRTGeometryHelper: public geo::AuxDetExptGeoHelperInterface {
  public:
    explicit CRTGeometryHelper(fhicl::ParameterSet const& pset);
  private:

    AuxDetChannelMapAlgPtr_t
    doConfigureAuxDetChannelMapAlg(fhicl::ParameterSet const & sortingParameters) const override;

    fhicl::ParameterSet fPset;
  };
}
DECLARE_ART_SERVICE_INTERFACE_IMPL(crt::CRTGeometryHelper, geo::AuxDetExptGeoHelperInterface, SHARED)

#endif // define
