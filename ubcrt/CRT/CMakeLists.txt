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
        CRTChannelMapAlg.cc 
        CRTGeometryHelper_service.cc
        CRTFileManager_service.cc
        CRTDetSim_module.cc
        CRTDaqSim_module.cc
        CRTHitProducer_module.cc
        CRTSimHitProducer_module.cc
        CRTSimHitCorr_module.cc
        CRTCalibration_module.cc
        CRTTrackProducer_module.cc
        CRTTzeroProducer_module.cc
	FlashExt_module.cc	
	CRTAuxFunctions.cc
  	CRTMerger_module.cc
	CRTFlashAna_module.cc
	TimeExtraction_module.cc
	TrackDump_module.cc
        CRTHitPositionPerModule_module.cc
	MergeTest_module.cc
	CRTRawtoCRTHit_source.cc
	CRTRawtoCRTHitGPScorr_source.cc
	CRTAnalyser_module.cc	
	CRTeffStd_module.cc			
	CRTTiming_module.cc
        T0recoCRT_module.cc                     
        T0recoCRTHit_module.cc                     
        T0recoCRTAnal_module.cc 
	CRTTiming_module.cc
        T0recoCRTHitAnal_module.cc 
    LIBRARY_NAME  
        CRTData
    LIB_LIBRARIES
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
	lardataobj_RecoBase
	lardataobj_AnalysisBase
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        cetlib cetlib_except
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
        ${MF_MESSAGELOGGER}
        
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        cetlib cetlib_except
)

simple_plugin( CRTGeometryHelper service
        larcorealg_Geometry
        CRT
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        
        ${ROOT_BASIC_LIB_LIST}
        cetlib cetlib_except
)

simple_plugin( CRTFileManager service
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
	gallery
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
	ifdh
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${ROOT_BASIC_LIB_LIST}
        cetlib cetlib_except
	${Boost_DATE_TIME_LIBRARY}
    BASENAME_ONLY
)

simple_plugin( CRTDetSim module
       	ubobj_CRT
	CRTData
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nurandom_RandomUtils_NuRandomService_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}        
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        cetlib cetlib_except
)

simple_plugin(CRTCalibration "module"
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
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}
        
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        cetlib cetlib_except
    BASENAME_ONLY
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )


simple_plugin(CRTSimHitProducer "module"
          larcorealg_Geometry
        larcore_Geometry_Geometry_service
          lardataalg_DetectorInfo
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
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
	${ART_ROOT_IO_TFILESERVICE_SERVICE}
        ${MF_MESSAGELOGGER}        
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        cetlib cetlib_except
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
        nutools_RandomUtils_NuRandomService_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}        
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        cetlib cetlib_except
    BASENAME_ONLY
)


simple_plugin(CRTHitProducer "module"
            CRT
              bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} 
	  cetlib cetlib_except
)

simple_plugin(CRTTrackProducer "module"
            CRT
              bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} 
	  cetlib cetlib_except
)
simple_plugin(CRTTzeroProducer "module"
            CRT
              bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} 
	  cetlib #cetlib_except
)

simple_plugin(FlashExt "module"
            CRT
              bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(TimeExtraction "module"
            CRT
	  bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(MergeTest "module"
            CRT
	  bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(TrackDump "module"
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
           art_root_io_tfile_support ${ROOT_CORE}
           art_Framework_Services_Registry
           art_root_io_TFileService_service
           ${ROOT_BASIC_LIB_LIST} )



simple_plugin(CRTRawtoCRTHit "source"
            CRT
	  bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTRawtoCRTHitGPScorr "source"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTAnalyser "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTeffStd "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(T0recoCRT "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(T0recoCRTHit "module"
            CRT
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(T0recoCRTAnal "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTTiming "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )
          
simple_plugin(CRTTriggerTimeProducer "module"
            CRT
          bernfebdaq-core_Overlays
          art_root_io_tfile_support ${ROOT_CORE}
          art_Framework_Services_Registry
          art_root_io_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(T0recoCRTHitAnal "module"
            CRT
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTTiming "module"
            CRT
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin( CRTMerger "module"
  CRT
  bernfebdaq-core_Overlays
  lardataobj_RecoBase
  CRTFileManager_service
  ${CETLIB}	
  ${ROOT_BASIC_LIB_LIST}
  ${ART_UTILITIES}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
  ${ARTDAQ_CORE}
  ${FHICLCPP}
  ${MF_MESSAGELOGGER}
  ${GALLERY}
  ubobj_CRT
  gallery
  cetlib cetlib_except	
)

simple_plugin( CRTFlashAna "module"
  CRT
  lardataobj_RecoBase
  ubobj_RawData
  ubobj_Trigger
  ubobj_CRT
  ${ART_UTILITIES}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${FHICLCPP}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
)


install_headers()
install_fhicl()
install_source()


