cet_make_library(
  SOURCE
  CRTAuxFunctions.cc
  CRTGeoObjectSorter.cc
  LIBRARIES
  PUBLIC
  larcorealg::Geometry
  messagefacility::MF_MessageLogger
  artdaq_core::artdaq-core_Data
)

cet_build_plugin(
  CRTAuxDetInitializer art::tool
  LIBRARIES
  PRIVATE
  larcorealg::Geometry
)

cet_build_plugin(
  CRTGeoObjectSorter art::tool
  LIBRARIES
  PRIVATE
  ubcrt::CRT
  larcorealg::Geometry
)

cet_build_plugin(
  CRTAnalyser art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(
  CRTCalibration art::EDAnalyzer
  LIBRARIES
  PRIVATE
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
  ROOT::Hist
  ROOT::Tree
)

#cet_build_plugin(
#  CRTDaqSim art::EDProducer
#  LIBRARIES
#  PRIVATE
#)

cet_build_plugin(
  CRTDataHitCorr art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  CRTDataQuality art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
  ROOT::Tree
)

cet_build_plugin(
  CRTDetSim art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  nurandom::RandomUtils_NuRandomService_service
)

cet_build_plugin(
  CRTFlashAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::Trigger
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Hist
)

cet_build_plugin(
  CRTHitPositionPerModule art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
)

cet_build_plugin(
  CRTHitProducer art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  ubcrt::CRT_CRTBernFEBDAQCore_Overlays
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CRTMerger art::EDProducer
  LIBRARIES
  PRIVATE
  ubcrt::CRT_CRTFileManager_service
  gallery::gallery
)

cet_build_plugin(
  CRTSimHitCorr art::EDProducer
  LIBRARIES
  PRIVATE
  ubevt::Database
  ubobj::CRT
  larevt::CalibrationDBI_Providers
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  art_root_io::TFileService_service
  nurandom::RandomUtils_NuRandomService_service
)

cet_build_plugin(
  CRTSimHitProducer art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::Simulation
  art_root_io::TFileService_service
)

cet_build_plugin(
  CRTTiming art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  art_root_io::TFileService_service
)

cet_build_plugin(
  CRTTrackProducer art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
)

#cet_build_plugin(
#  CRTTriggerFilter
#  LIBRARIES
#  PRIVATE
#)

cet_build_plugin(
  CRTTriggerTimeProducerSimple art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  larpandora::LArPandoraInterface
  larcore::ServiceUtil
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CRTTriggerTimeProducer art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  larpandora::LArPandoraInterface
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CRTTzeroProducer art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
)

cet_build_plugin(
  CRTeffStd art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
)

cet_build_plugin(
  FlashExt art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
  ROOT::Tree
)

cet_build_plugin(
  MergeTest art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
)

cet_build_plugin(
  T0recoCRTAnal art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  lardata::Utilities
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  T0recoCRTHitAna2 art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::SpaceCharge
  lardata::DetectorClocksService
  larcore::ServiceUtil
  lardataobj::RecoBase
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
)

cet_build_plugin(
  T0recoCRTHitAnal art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  T0recoCRTHit art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::SpaceCharge
  lardata::DetectorPropertiesService
  larcore::ServiceUtil
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
)

cet_build_plugin(
  T0recoCRTTrackAnal art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  lardata::Utilities
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  T0recoCRTTrack art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
)

#cet_build_plugin(
#  T0recoCRT art::EDProducer
#  LIBRARIES
#  PRIVATE
#)

cet_build_plugin(
  TimeExtraction art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  artdaq_core::artdaq-core_Data
  art_root_io::TFileService_service
)

cet_build_plugin(
  TrackDump art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::CRT
  larevt::SpaceCharge
  lardata::Utilities
  larcore::ServiceUtil
  lardataobj::RecoBase
  art_root_io::TFileService_service
  artdaq_core::artdaq-core_Data
  ROOT::Tree
)

cet_build_plugin(
  CRTFileManager art::service
  LIBRARIES
  PUBLIC
  ubobj::CRT
  ifdh_art::IFDH_service
  gallery::gallery
)

cet_build_plugin(
  CRTRawtoCRTHit art::source
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CRTRawtoCRTHitGPScorr art::source
  LIBRARIES
  PRIVATE
  ubobj::CRT
  ubcrt::CRT
  art_root_io::TFileService_service
  ROOT::Tree
)

add_subdirectory(CRTBernFEBDAQCore)
#add_subdirectory(CRTBinarySource)

install_headers()
install_fhicl()
install_source()
