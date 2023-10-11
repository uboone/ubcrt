cet_build_plugin(
  CRTDistanceFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  ExternalVtx art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  G4VetoSignalImpact art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larpandora::LArPandoraInterface
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  UBCRTCosmicFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

art_dictionary(
  DICTIONARY_LIBRARIES
  ubobj::CRT
  lardata::Utilities
)

install_headers()
install_fhicl()
install_source()
