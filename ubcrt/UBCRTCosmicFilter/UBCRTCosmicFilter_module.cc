////////////////////////////////////////////////////////////////////////
// Class:       UBCRTCosmicFilter
// Plugin Type: filter (art v2_11_03)
// File:        UBCRTCosmicFilter_module.cc
//
// Generated at Fri Sep 28 13:12:47 2018 by Christopher Barnes using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// declare the package with the CRT hit information.
#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

#include <memory>

class UBCRTCosmicFilter;


class UBCRTCosmicFilter : public art::EDFilter {
public:
  explicit UBCRTCosmicFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBCRTCosmicFilter(UBCRTCosmicFilter const &) = delete;
  UBCRTCosmicFilter(UBCRTCosmicFilter &&) = delete;
  UBCRTCosmicFilter & operator = (UBCRTCosmicFilter const &) = delete;
  UBCRTCosmicFilter & operator = (UBCRTCosmicFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:

  // Declare member data here.
  std::string fBeamFlashProducer;
  std::string fCRTHitProducer;
  std::string fDAQHeaderProducer;
  double      fBeamStart;
  double      fBeamEnd;
  double      fPEMin;
  double      fDTOffset;
  double      fResolution;
  bool        fSingleBeamgateFlashOnly;
  bool        verbose;

  // beam info
  double      beam_flash_time;
  double      beam_flash_PEs;
  
  // CRT hit time.
  double      closest_CRT_hit_time;
  
  // Counter for the events.
  int    event_counter;

};


UBCRTCosmicFilter::UBCRTCosmicFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  fBeamFlashProducer       = "simpleFlashBeam";
  fCRTHitProducer          = "merger";
  fDAQHeaderProducer       = "daq";
  fBeamStart               = 5.0;
  fBeamEnd                 = 16.0;
  fPEMin                   = 50.0;
  fDTOffset                = 68600.0;
  fResolution              = 1.0;
  fSingleBeamgateFlashOnly = true;
  verbose                  = true;
}

bool UBCRTCosmicFilter::filter(art::Event & e)
{
  // Implementation of required member function here.

  std::cout << "Now looping over event #" << event_counter << "." << std::endl;
  event_counter++;

  // Declare an object for the 
  art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;
  e.getByLabel(fDAQHeaderProducer, rawHandle_DAQHeader);

  raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
  double evt_timeGPS_nsec = evtTimeGPS.timeLow();

  // load CRT hits.                                                                                                                                                                                    
  art::Handle<std::vector<crt::CRTHit> > crthit_h;
  e.getByLabel( fCRTHitProducer , crthit_h );

  // make sure CRT hits look good.                                                                                                                                                                   
  if(!crthit_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate CRT Hit!"<<std::endl;
    throw std::exception();
  }

  // load beam flashes here.                                                                                                                                                                           
  art::Handle<std::vector<recob::OpFlash> > beamflash_h;
  e.getByLabel(fBeamFlashProducer , beamflash_h);

  // make sure beam flashes look good.                                                                                                                                                                  
  if(!beamflash_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Beam Flash!"<<std::endl;
    throw std::exception();
  }

  // Loop through the beam flashes and find the closest crt hit to the beam flashes.                                                                                                                    
  if ( verbose ) 
    std::cout << "Number of beam flashes in this event = " << beamflash_h->size() << "." << std::endl;

  double closest_CRT_diff     = 100000.0;
  closest_CRT_hit_time        = 0.;

  bool within_resolution = false;

  for ( size_t i = 0; i < beamflash_h->size(); i++ ) {

    beam_flash_time = beamflash_h->at( i ).Time();
    beam_flash_PEs  = beamflash_h->at( i ).TotalPE();

    if ( verbose ) {
      std::cout << "Time of the beam flash = " << beam_flash_time << " us." << std::endl;
      std::cout << "PEs of the beam flash = " << beam_flash_PEs << " PEs." << std::endl;
    }

    // Only look at those flashes that are in the beamgate and of the minimum amount of PEs.
    if ( beam_flash_time < fBeamStart || beam_flash_time > fBeamEnd || beam_flash_PEs < fPEMin ) {
      
      if ( verbose ) 
	std::cout << "Flash is outside of beam gate or too low in intensity.  Skipping!" << std::endl;
      
      continue;
    
    }

    // Loop over the CRT hits.                                                                                                                                                                           
    for ( size_t j = 0; j < crthit_h->size(); j++ ) {

      if ( verbose ) 
	std::cout << "Time of the CRT Hit wrt the event timestamp = " << ( ( crthit_h->at( j ).ts0_ns - evt_timeGPS_nsec + fDTOffset ) / 1000. ) << " us." << std::endl;

      // Reset the variables.                                                                                                                                                                           
      if ( fabs( beam_flash_time - ( ( crthit_h->at( j ).ts0_ns - evt_timeGPS_nsec  + fDTOffset ) / 1000. ) ) < closest_CRT_diff ) {
        closest_CRT_diff     = fabs( beam_flash_time - ( ( crthit_h->at( j ).ts0_ns - evt_timeGPS_nsec + fDTOffset ) / 1000. ) );
        closest_CRT_hit_time = ( ( crthit_h->at( j ).ts0_ns - evt_timeGPS_nsec +  fDTOffset ) / 1000. );

	// set 'within_resolution' to 'true' and break the loop if 'closest_crt_diff' is less than fResolution.
	if ( closest_CRT_diff < fResolution ) {
	  within_resolution = true;
	  break;
	}

      } // End of conditional for closest CRT hit time.                                                                                                                                                   

    } // End of loop over CRT hits.                                                                                                                                                                         
    // End the loop over beam flashes if 'within_resolution' = 'true'.  You do not have to consider any more flashes if that is the case.
    if ( within_resolution == true ) 
      break;

  } // End of loop over beam flashes.                                                                                                                                                                       
  
  // Look at cases in which you could return true.
  if ( within_resolution == true ) {

    // If you want a single flash and there is more than one in the beamgate, than return false.
    if ( fSingleBeamgateFlashOnly == true and beamflash_h->size() != 1 ) 
      return false;

    // In all other cases return true.
    else
      return true;

  }

  // If it is not within resolution, return false.
  else 
    return false;

} // End of the filter module
DEFINE_ART_MODULE(UBCRTCosmicFilter)
