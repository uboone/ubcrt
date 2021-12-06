art_make_library ( LIBRARY_NAME bernfebdaq-core_Overlays
		   SOURCE BernFEBFragment.cc BernZMQFragment.cc
		   LIBRARIES
		   artdaq_core::artdaq-core_Utilities
  	    	   cetlib::cetlib
  	    	   cetlib_except::cetlib_except
  		 )
install_headers()
install_source()
