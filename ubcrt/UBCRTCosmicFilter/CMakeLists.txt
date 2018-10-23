if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-return-type-c-linkage )
endif()

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY

          LIB_LIBRARIES
                        larpandora_LArPandoraInterface
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        
          MODULE_LIBRARIES 
                        larpandora_LArPandoraInterface
                        lardataobj_RecoBase	
                        lardataobj_RecoBase
                        nusimdata_SimulationBase		
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}                        
                        ${FHICLCPP}
                        cetlib cetlib_except
			nusimdata_SimulationBase
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}

        )

install_headers()
install_fhicl()
install_source()
