cet_make_library(
  SOURCE
  BernZMQBinaryInputStreamReader.cc
  LIBRARIES
  PUBLIC
  ubcrt::CRT_CRTBernFEBDAQCore_Overlays
  canvas::canvas
)

cet_build_plugin(
  BernZMQBinaryInput art::source
  LIBRARIES
  PRIVATE
  ubcrt::CRT_CRTBernFEBDAQCore_BernZMQBinarySource
)

install_headers()
install_source()
