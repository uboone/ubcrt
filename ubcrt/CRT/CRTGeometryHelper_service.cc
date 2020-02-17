#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeo.h"

#include "ubcrt/CRT/CRTGeometryHelper.hh"
#include "ubcrt/CRT/CRTChannelMapAlg.hh"

#include <vector>

namespace crt
{

  CRTGeometryHelper::CRTGeometryHelper( fhicl::ParameterSet const & pset ) : fPset( pset ){}

  CRTGeometryHelper::AuxDetChannelMapAlgPtr_t
  CRTGeometryHelper::doConfigureAuxDetChannelMapAlg(fhicl::ParameterSet const & sortingParameters) const
  {
    return std::make_unique<CRTChannelMapAlg>( fPset, sortingParameters );
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(crt::CRTGeometryHelper, geo::AuxDetExptGeoHelperInterface)
