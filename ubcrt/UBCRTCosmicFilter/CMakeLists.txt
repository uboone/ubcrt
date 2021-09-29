if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-return-type-c-linkage )
endif()

art_make( BASENAME_ONLY
          DICT_LIBRARIES
                        ${ART_UTILITIES}
                        messagefacility::MF_MessageLogger
          MODULE_LIBRARIES
                        larpandora_LArPandoraInterface
                        lardataobj_RecoBase
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        ${MF_MESSAGELOGGER}
                        nusimdata_SimulationBase
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}


        )

install_headers()
install_fhicl()
install_source()
