cet_make_library(
  SOURCE
  BernFEBFragment.cc
  BernZMQFragment.cc
  LIBRARIES
  PUBLIC
  artdaq_core::artdaq-core_Data
)

install_headers()
install_source()
