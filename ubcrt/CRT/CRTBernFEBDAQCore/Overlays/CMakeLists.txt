art_make_library ( LIBRARY_NAME bernfebdaq-core_Overlays
		   SOURCE BernFEBFragment.cc BernZMQFragment.cc
		   LIBRARIES
  	    	   ${ARTDAQ_DAQDATA}
		   ${ARTDAQ-CORE_UTILITIES}
  	    	   ${CETLIB}
  	    	   ${CETLIB_EXCEPT}
  		 )
install_headers()
install_source()
