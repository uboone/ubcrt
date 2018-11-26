if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-dynamic-exception-spec -Wno-return-type-c-linkage -Wno-undefined-var-template )
endif()

# canvas_utlities
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )

art_make( BASENAME_ONLY
    EXCLUDE 
        CRTAnaFun.cc
	CRTXSEC_module.cc 
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
        CRTAna
    SOURCE
        CRTAnaFun.cc 
    LIBRARIES 
        #larcorealg_Geometry
        ${MF_MESSAGELOGGER}
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        cetlib cetlib_except
)


simple_plugin(CRTXSEC "module"
            CRTAna
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )


install_headers()
install_fhicl()
install_source()


