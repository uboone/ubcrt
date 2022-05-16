if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-return-type-c-linkage )
endif()

art_make( BASENAME_ONLY
          DICT_LIBRARIES
                        art::Utilities
                        messagefacility::MF_MessageLogger
          MODULE_LIBRARIES
                        larpandora_LArPandoraInterface
                        lardataobj_RecoBase
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        messagefacility::MF_MessageLogger
                        nusimdata::SimulationBase
                        CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
                        ROOT::Tree
        )

install_headers()
install_fhicl()
install_source()
