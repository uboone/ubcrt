////////////////////////////////////////////////////////////////////////
// Class:       UBProximityFilter
// Plugin Type: filter (art v2_11_03)
// File:        UBProximityFilter_module.cc
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
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// declare the package with the CRT hit information.
#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// Declare the file for the associations to be made.
#include "lardata/Utilities/AssociationUtil.h"

// C++
#include <memory>

// ROOT
#include <TTree.h>

class UBProximityFilter;

class UBProximityFilter : public art::EDFilter
{
public:
  explicit UBProximityFilter(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBProximityFilter(UBProximityFilter const &) = delete;
  UBProximityFilter(UBProximityFilter &&) = delete;
  UBProximityFilter &operator=(UBProximityFilter const &) = delete;
  UBProximityFilter &operator=(UBProximityFilter &&) = delete;

  // Use the 'beginJob()' function to save information in the TTree.
  void beginJob() override;

  // Required functions.
  bool filter(art::Event &e) override;

  // Optional functions.
  float calculateDistance (float x0, float y0, float z0, float x1, float y1, float z1);
private:
  // Declare member data here.
  std::string fTrackProducer;
  std::string fCRTHitProducer;
  float fDistanceSelectionValue;
  bool  fCheckCRTHitAssn;
  bool  fCheckACPTAssn;
  bool  fuseAsFilter;
  bool  verbose;
 

  // TTree Declaration.
  TTree *_tree;

  // Event info.
  int   _run;
  int   _subrun;
  int   _event;
  float _minTrkVtxDistance;
  bool  _keepEvent;
  float vtx_X; float vtx_Y; float vtx_Z;

  // Counter for the events.
  int event_counter;
};

UBProximityFilter::UBProximityFilter(fhicl::ParameterSet const &p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  fTrackProducer  = p.get<std::string>("TrackProducer");
  fCRTHitProducer = p.get<std::string>("CRTHitProducer");
  fuseAsFilter    = p.get<bool>("useAsFilter");
  verbose         = p.get<bool>("verbose");
  fDistanceSelectionValue = p.get<std::string>("DistanceSelectionValue");
  fCheckCRTHitAssn        = p.get<bool>("CheckCRTHitAssn");
  fCheckACPTAssn          = p.get<bool>("CheckACPTAssn");
}

void UBProximityFilter::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree", "Flash/CRT Matching Information For All Events");
  _tree->Branch("_run"   , &_run   , "run/I");
  _tree->Branch("_subrun", &_subrun, "subrun/I");
  _tree->Branch("_event" , &_event , "event/I");
}


bool UBProximityFilter::calculateDistance(float x0, float y0, float z0, float x1, float y1, float z1)
{
  return TMath::Sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1) );
}
bool UBProximityFilter::filter(art::Event &e)
{
  bool withinDistance = false;
  // Implementation of required member function here.

  // Fill the event info.
  _event  = e.event();
  _subrun = e.subRun();
  _run    = e.run();

  if (verbose)
    std::cout << "Now looping over event #" << event_counter << "." << std::endl;

  // Iterate the event counter.
  event_counter++;

  // Take the neutrino vtx
  

  // Take the tracks
  art::Handle< std::vector<recob::Track>  > trackListHandle; 
  std::vector<art::Ptr<recob::Track> >  tracklist;
  if (!evt.getByLabel(fTrackProducer,trackListHandle)) continue;
  art::fill_ptr_vector(tracklist, trackListHandle);
  //check to make sure the data we asked for is valid
  if(!trackListHandle.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
	      << ", event " << evt.event() << " has zero"
	      << " tracks " << " in module " << fTrackModuleLabel << std::endl;
    std::cout << std::endl;
    return;
  }


  
  // Loop on the tracklist
  for(size_t j = 0; j < tracklist.size();; j++) {    
    art::Ptr<recob::Track> ptrack(trackListHandle, j);
    const recob::Track& track = *ptrack;

    
    for (size_t iPt = 0; iPt < track.NumberTrajectoryPoints(); iPt++ )
      {
	auto thisTrajPointPosition = track.TrajectoryPoint (iPt).position;

	float thisPointDist = calculateDistance(thisTrajPointPosition.X(), thisTrajPointPosition.Y(), thisTrajPointPosition.Z(), 
						nuvtx_X, nuvtx_Y, nuvtx_Z);

      }


    if (fCheckCRTHitAssn)
      {
      }
    if (fCheckACPTAssn)
      {}
    //art::FindMany<anab::T0> trk_t0A_assn_v(trackListHandle, evt, data_label_t0A_);
    //const std::vector<const anab::T0*>& T0_acpt = trk_t0A_assn_v.at(j);
    
  }




  // Take the neutrino vtx
  // Calculate impact parameter
 
 
  // Make the cut
  if (_minTrkVtxDistance < fDistanceSelectionValue) withinDistance = true;
  
  _keepEvent == false;
  _tree->Fill();




  // Return true if you are within resolution.
  if (fuseAsFilter && _keepEvent == false)
    return false;

  // Otherwise return false.
  return true;

} // End of the filter module


DEFINE_ART_MODULE(UBProxityFilter)
