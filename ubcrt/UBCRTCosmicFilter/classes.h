// Structural Includes.
#include "canvas/Persistency/Common/Wrapper.h"
#include "lardata/Utilities/AssociationUtil.h"

// Data Includes.
#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"

template class art::Assns<crt::CRTHit,recob::OpFlash,void>;
template class art::Assns<recob::OpFlash,crt::CRTHit,void>;
template class art::Wrapper<art::Assns<crt::CRTHit,recob::OpFlash,void> >;
template class art::Wrapper<art::Assns<recob::OpFlash,crt::CRTHit,void> >;
