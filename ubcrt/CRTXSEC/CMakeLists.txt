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
        CRTEventAnalyser_module.cc
        CRTtoTTree_module.cc
        CRTtoTTreeVector_module.cc
        EventAnalyser_module.cc
        CRTTPCAssoAnalyser_module.cc
        TPCTrackAnalyser_module.cc
        CRTLArPandoraExternalEventBuilding_module.cc
        CRTSimpleNeutrinoId_tool.cc
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
  lardata_Utilities
  larpandora_LarPandoraEventBuilding
  larpandora_LarPandoraInterface
  larpandora_LArPandoraObjects
  larpandoraobj
  LArPandoraContent
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}   
        cetlib cetlib_except
)

art_make_library( 
    LIBRARY_NAME 
        CRTAna
    SOURCE
        CRTAnaFun.cc      
    LIBRARIES 
        #larcorealg_Geometry
        #CRTData
        CRT
        bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          lardataobj_RecoBase
          ##lardataobj_AnalysisBase
          ${LARDATA_LIB}
          lardata_Utilities
        ${MF_MESSAGELOGGER}
        ${ART_FRAMEWORK_CORE}
		    ${ART_FRAMEWORK_PRINCIPAL}
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}   
        cetlib cetlib_except
)


simple_plugin(CRTXSEC "module"
            CRT
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )
          
simple_plugin(CRTEventAnalyser "module"
            CRT
            CRTAna
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )
          
simple_plugin(CRTtoTTree "module"
            CRT
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )
          
simple_plugin(CRTtoTTreeVector "module"
            CRT
            CRTAna
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )
#[[          
##simple_plugin(EventAnalyser "module"
            CRT
#            CRTAna
#          bernfebdaq-core_Overlays
#          LArPandoraContent
#          larpandoraobj
#          art_Framework_Services_Optional
#          art_Framework_Services_Registry
#          art_Framework_Services_Optional_TFileService_service
#          ${ROOT_BASIC_LIB_LIST} )
]]        
simple_plugin(CRTTPCAssoAnalyser "module"
            CRT
            CRTAna
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )
          
simple_plugin(TPCTrackAnalyser "module"
            CRT
            CRTAna
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )
 #[[         
simple_plugin(CRTSimpleNeutrinoId "tool" larpandora_LArPandoraEventBuilding)
          
simple_plugin(CRTLArPandoraExternalEventBuilding "module"
            CRT
            CRTAna
          bernfebdaq-core_Overlays
          lardataobj_RecoBase
          lardataobj_AnalysisBase
          lardata_Utilities
          #larpandora_LarPandoraEventBuilding
          #larpandora_LarPandoraInterface
          #larpandora_LArPandoraObjects
          #larpandoraobj
          #LArPandoraContent
         art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )
]]
install_headers()
install_fhicl()
install_source()


