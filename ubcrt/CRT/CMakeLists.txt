if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-dynamic-exception-spec -Wno-return-type-c-linkage -Wno-undefined-var-template )
endif()

# canvas_utlities
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )

add_subdirectory(CRTBernFEBDAQCore)
#add_subdirectory(CRTBinarySource)

art_make( BASENAME_ONLY
    EXCLUDE
        CRTTriggerTimeProducer_module.cc
        CRTTriggerTimeProducerSimple_module.cc
        CRTChannelMapAlg.cc 
        CRTGeometryHelper_service.cc
        CRTFileManager_service.cc
        CRTDetSim_module.cc
        CRTDaqSim_module.cc
        CRTHitProducer_module.cc
        CRTSimHitProducer_module.cc
        CRTSimHitCorr_module.cc
        CRTDataHitCorr_module.cc
        CRTCalibration_module.cc
        CRTTrackProducer_module.cc
        CRTTzeroProducer_module.cc
	FlashExt_module.cc	
	CRTAuxFunctions.cc
  	CRTMerger_module.cc
	CRTFlashAna_module.cc
	TimeExtraction_module.cc
	TrackDump_module.cc
	CRTDataQuality_module.cc
        CRTHitPositionPerModule_module.cc
	MergeTest_module.cc
	CRTRawtoCRTHit_source.cc
	CRTRawtoCRTHitGPScorr_source.cc
	CRTAnalyser_module.cc	
	CRTeffStd_module.cc			
	CRTTiming_module.cc
        T0recoCRT_module.cc
	T0recoCRTTrack_module.cc                     
        T0recoCRTHit_module.cc                     
        T0recoCRTAnal_module.cc 
	T0recoCRTTrackAnal_module.cc
        T0recoCRTHitAnal_module.cc 
	T0recoCRTHitAna2_module.cc 
    LIBRARY_NAME  
        CRTData
    LIB_LIBRARIES
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Persistency_Common canvas::canvas
        art::Persistency_Provenance canvas::canvas
        art::Utilities canvas::canvas
        lardataobj_RecoBase
        lardataobj_AnalysisBase
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        fhiclcpp::fhiclcpp
        ROOT::Geom
        cetlib::cetlib cetlib_except
)

art_make_library(
    LIBRARY_NAME
        CRT
    SOURCE
        CRTChannelMapAlg.cc
        CRTAuxFunctions.cc
    LIBRARIES
        #larcorealg_Geometry
        CRTData
        messagefacility::MF_MessageLogger

        fhiclcpp::fhiclcpp

        ROOT::Geom
        cetlib::cetlib cetlib_except
)

simple_plugin( CRTGeometryHelper service
        larcorealg_Geometry
        CRT
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Persistency_Common canvas::canvas
        art::Persistency_Provenance canvas::canvas
        art::Utilities canvas::canvas
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger


        cetlib::cetlib cetlib_except
)

simple_plugin( CRTFileManager service
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        gallery
        art::Persistency_Common canvas::canvas
        art::Persistency_Provenance canvas::canvas
        art::Utilities canvas::canvas
        ifdh::ifdh
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        

        cetlib::cetlib cetlib_except
        Boost::date_time
    BASENAME_ONLY
)

simple_plugin( CRTDetSim module
        ubobj_CRT
        CRTData
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Persistency_Common canvas::canvas
        art::Persistency_Provenance canvas::canvas
        art::Utilities canvas::canvas
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger

        CLHEP::CLHEP
        cetlib::cetlib cetlib_except
)

simple_plugin(CRTCalibration "module"
        CRT
        CRTData
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Persistency_Common canvas::canvas
        art::Persistency_Provenance canvas::canvas
        art::Utilities canvas::canvas
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
      nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger


        CLHEP::CLHEP
        cetlib::cetlib cetlib_except
    BASENAME_ONLY
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
           )


simple_plugin(CRTSimHitProducer "module"
          larcorealg_Geometry
        larcore_Geometry_Geometry_service
          lardataalg_DetectorInfo
          CRT
          CRTData
                  art::Framework_Core
                  art::Framework_IO_Sources
                  art::Framework_Principal
        art::Persistency_Common canvas::canvas
        art::Persistency_Provenance canvas::canvas
        art::Utilities canvas::canvas
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        art_root_io::TFileService_service
        messagefacility::MF_MessageLogger

        CLHEP::CLHEP
        cetlib::cetlib cetlib_except
)

simple_plugin(CRTSimHitCorr "module"
          larcorealg_Geometry
        larcore_Geometry_Geometry_service
	  CRT
	  CRTData
	          ${ART_FRAMEWORK_CORE}
        	  ${ART_FRAMEWORK_IO_SOURCES}
        	  ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom_RandomUtils_NuRandomService_service
        art_root_io_tfile_support ${ROOT_CORE}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        art_root_io_TFileService_service
        ${MF_MESSAGELOGGER}        
        ${ROOT_BASIC_LIB_LIST}
        CLHEP::CLHEP
        cetlib cetlib_except
    BASENAME_ONLY
)


simple_plugin(CRTDataHitCorr "module"
          larcorealg_Geometry
        larcore_Geometry_Geometry_service
	  CRT
	  CRTData
	          ${ART_FRAMEWORK_CORE}
        	  ${ART_FRAMEWORK_IO_SOURCES}
        	  ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom_RandomUtils_NuRandomService_service
 	art_root_io_tfile_support ${ROOT_CORE}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        art_root_io_TFileService_service
        ${MF_MESSAGELOGGER}        
        ${ROOT_BASIC_LIB_LIST}
        CLHEP::CLHEP
        cetlib cetlib_except
    BASENAME_ONLY
)

simple_plugin(CRTHitProducer "module"
            CRT
              bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          ROOT::Tree
          cetlib::cetlib
)

simple_plugin(CRTTrackProducer "module"
            CRT
              bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          cetlib::cetlib
)
simple_plugin(CRTTzeroProducer "module"
            CRT
              bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          cetlib::cetlib
)

simple_plugin(FlashExt "module"
            CRT
              bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          ROOT::Tree
           )

simple_plugin(TimeExtraction "module"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
           )

simple_plugin(MergeTest "module"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
           )

simple_plugin(TrackDump "module"
            CRT
          bernfebdaq-core_Overlays
          lardataobj_RecoBase
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          ROOT::Tree
           )

simple_plugin(CRTDataQuality "module"
            CRT
	  bernfebdaq-core_Overlays
	  lardataobj_RecoBase
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
	  art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTHitPositionPerModule "module"
             CRT
           bernfebdaq-core_Overlays
           lardataobj_RecoBase
           art::Framework_Services_Registry
           art_root_io::TFileService_service
            )



simple_plugin(CRTRawtoCRTHit "source"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          ROOT::Tree
           )

simple_plugin(CRTRawtoCRTHitGPScorr "source"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          ROOT::Tree
           )

simple_plugin(CRTAnalyser "module"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          ROOT::Tree
           )

simple_plugin(CRTeffStd "module"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          nusimdata::SimulationBase
           )

simple_plugin(T0recoCRT "module"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          nusimdata::SimulationBase
           )

simple_plugin(T0recoCRTTrack "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
	  art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(T0recoCRTTrackAnal "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
	  art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(T0recoCRTHit "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
	  art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(T0recoCRTAnal "module"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          nusimdata::SimulationBase
           )

simple_plugin(CRTTiming "module"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
           )

simple_plugin(CRTTriggerTimeProducer "module"
            CRT
          bernfebdaq-core_Overlays
          art::Framework_Services_Registry
          art_root_io::TFileService_service
          ROOT::Tree
           )

simple_plugin(T0recoCRTHitAnal "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
	  art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )


simple_plugin(T0recoCRTHitAna2 "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
	  art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTTriggerTimeProducerSimple "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
	  art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )


simple_plugin( CRTMerger "module"
  CRT
  bernfebdaq-core_Overlays
  lardataobj_RecoBase
  CRTFileManager_service
  cetlib::cetlib

  art::Utilities
  art::Framework_Services_Registry
  ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
  ##${ARTDAQ_CORE}
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  ${GALLERY}
  ubobj_CRT
  gallery
  cetlib::cetlib cetlib_except
)

simple_plugin( CRTFlashAna "module"
  CRT
  lardataobj_RecoBase
  ubobj_RawData
  ubobj_Trigger
  ubobj_CRT
  art::Utilities
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  fhiclcpp::fhiclcpp
  art::Persistency_Common canvas::canvas
  art::Persistency_Provenance canvas::canvas
  art::Utilities canvas::canvas
  messagefacility::MF_MessageLogger
  
  cetlib::cetlib cetlib_except

)


install_headers()
install_fhicl()
install_source()
