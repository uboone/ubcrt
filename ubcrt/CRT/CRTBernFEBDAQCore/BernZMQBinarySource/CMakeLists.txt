art_make( LIBRARY_NAME bernfebdaq_BernZMQBinarySource
	  	 LIB_LIBRARIES
				bernfebdaq-core_Overlays
				fhiclcpp::fhiclcpp
				art::Persistency_Provenance				
	  	 SOURCE_LIBRARIES
				bernfebdaq_BernZMQBinarySource
		 		bernfebdaq-core_Overlays
				art::Framework_Core
				art::Framework_IO_Sources
				art::Framework_Principal
				art::Persistency_Provenance
				canvas::canvas
				${CANVAS_LIB_LIST}
				messagefacility::MF_MessageLogger
				${ROOT_CINTEX}
				ROOT::RIO
				ROOT::Core
				ROOT::Net
				
				cetlib::cetlib
				cetlib_except::cetlib_except
				${BOOST_LIB}
				fhiclcpp::fhiclcpp
	         BASENAME_ONLY
	       )

install_headers()
install_source()
