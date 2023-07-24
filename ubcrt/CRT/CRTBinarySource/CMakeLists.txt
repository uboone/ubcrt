art_make( LIBRARY_NAME crt_CRTRawBinarySource
                 LIB_LIBRARIES
                                fhiclcpp::fhiclcpp
                                art::Persistency_Provenance                              
                 SOURCE_LIBRARIES
                                art::Framework_Core
                                art::Framework_IO_Sources
                                art::Framework_Principal
                                art::Persistency_Provenance
                                canvas::canvas
                                messagefacility::MF_MessageLogger
                                ROOT::RIO
                                ROOT::Core
                                ROOT::Net
                                
                                cetlib::cetlib
                                cetlib_except::cetlib_except
                                fhiclcpp::fhiclcpp
                 BASENAME_ONLY
               ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.21.01) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types



install_headers()
install_source()
