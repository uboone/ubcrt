////////////////////////////////////////////////////////////////////////
// Class:       CRTDistanceFilter
// File:        CRTDistanceFilter_module.cc
//
// Generated at Wed June 26 2019 by Elena Gramellini elenag@fnal.gov
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
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"


// include track information
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
//#include "lardataobj/RecoBase/TrackHitMeta.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// Declare the file for the associations to be made.
#include "lardata/Utilities/AssociationUtil.h"

// C++
#include <memory>
#include <array>

// ROOT
#include <TTree.h>

class CRTDistanceFilter;

class CRTDistanceFilter : public art::EDFilter
{
public:
  explicit CRTDistanceFilter(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTDistanceFilter(CRTDistanceFilter const &) = delete;
  CRTDistanceFilter(CRTDistanceFilter &&) = delete;
  CRTDistanceFilter &operator=(CRTDistanceFilter const &) = delete;
  CRTDistanceFilter &operator=(CRTDistanceFilter &&) = delete;

  // Use the 'beginJob()' function to save information in the TTree.
  void beginJob() override;
  void ResetVar();
  // Required functions.
  bool filter(art::Event &e) override;
  void readExtVertices( std::vector<art::Ptr<recob::Vertex> > & vtxlist );

private:
  // Declare member data here.
  std::string fTrackModuleLabel;
  std::string fVtxModuleLabel;
  std::string fCRTHitLabel;
  std::string fCRTTrackAssnProducer;


  bool  fuseAsFilter;
  bool  fVerbose;
  bool  fIsThisMC;
  bool  fUseCRTTag;
  bool  fUseExternalVtx;
  float fMinDistance;
  /// 

  // TTree Declaration.
  TTree *_tree;

  // Event info.
  int   _run;
  int   _subrun;
  int   _event;
  int   _nVtxs;
  int   _nTracks;
  int   _nTrackCRTAssn;
  int   _nGoodVtx;
  bool  _keepEvent = false;
  std::vector<float> _closestDistanceVtxTrk;

};




CRTDistanceFilter::CRTDistanceFilter(fhicl::ParameterSet const &p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  produces< art::Assns<recob::Vertex, recob::Track> >();
  produces< std::vector<recob::Vertex>   >(); 

  fTrackModuleLabel     = p.get< std::string >("TrackModuleLabel"   , "pandora");
  fVtxModuleLabel       = p.get< std::string >("VtxModuleLabel"     , "pandora");
  //OverlayRecoStage2... | crttrackmatch............. | ..................... | art::Assns<recob::Track,crt::CRTHit,void>............................... | ...13
  fCRTHitLabel          = p.get<std::string>("CRTHitLabel"         ,"crthitcorr");
  fCRTTrackAssnProducer = p.get<std::string>("CRTTrackAssnProducer","crttrackmatch");
  fuseAsFilter          = p.get<bool >("useAsFilter");
  fVerbose              = p.get<bool >("verbose"       ,false);
  fIsThisMC             = p.get<bool >("isThisMC"      ,false);
  fUseCRTTag            = p.get<bool >("useCRTTag"     ,true);
  fUseExternalVtx       = p.get<bool >("useExternalVtx",false);
  fMinDistance          = p.get<float>("minDistance",14.);
}

void CRTDistanceFilter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree", "Track/Vtx Matching Information For All Events");
  _tree->Branch("_run"          , &_run             , "run/I");
  _tree->Branch("_subrun"       , &_subrun          , "subrun/I");
  _tree->Branch("_event"        , &_event           , "event/I");
  _tree->Branch("_nVtxs"        , &_nVtxs           , "nVtxs/I");
  _tree->Branch("_nTracks"      , &_nTracks         , "nTracks/I");
  _tree->Branch("_nTrackCRTAssn", &_nTrackCRTAssn   , "nTrackCRTAssn/I");
  _tree->Branch("_nGoodVtx"     , &_nGoodVtx        , "nGoodVtx/I");
  _tree->Branch("_keepEvent"    , &_keepEvent       , "keepEvent/O");
  _tree->Branch("_closestDistanceVtxTrk", "std::vector< float >", &_closestDistanceVtxTrk);   
}


void CRTDistanceFilter::ResetVar()
{
  _run              = -999999;
  _subrun           = -999999;
  _event            = -999999;
  _nVtxs            = -999999;
  _nTracks          = -999999;
  _nTrackCRTAssn    = -999999;
  _nGoodVtx         = -999999;
  _keepEvent        = -999999;
  _closestDistanceVtxTrk.clear();   


}


bool CRTDistanceFilter::filter(art::Event &e)
{
  // Clear the tree
  ResetVar();

  // This is what we want to produce
  std::unique_ptr<std::vector<recob::Vertex> > VtxCol(new std::vector<recob::Vertex>);
  std::unique_ptr<art::Assns<recob::Vertex, recob::Track>> vertex_track_assn_v(new art::Assns<recob::Vertex, recob::Track>);

  // Fill the event info.
  _event = e.event();
  _subrun = e.subRun();
  _run = e.run();

  // Get the reconstructed verteces
  art::Handle< std::vector<recob::Vertex> > vtxListHandle;
  std::vector<art::Ptr<recob::Vertex> > vtxlist;
  if (!e.getByLabel(fVtxModuleLabel,vtxListHandle)) {  e.put(std::move(VtxCol));  e.put(std::move(vertex_track_assn_v));  return false;}
  art::fill_ptr_vector(vtxlist, vtxListHandle);
    
  // Get the reconstructed tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (!e.getByLabel(fTrackModuleLabel,trackListHandle)) {  e.put(std::move(VtxCol)); e.put(std::move(vertex_track_assn_v)); return false;}
  art::fill_ptr_vector(tracklist, trackListHandle);


  // Are we going to use the CRT information?
  // if so, let's fetch which reco tracks have a CRT Hit association
  std::vector<int> matchedRecoTrkKey;
  matchedRecoTrkKey.clear();
  if (fUseCRTTag){
    // Get CRT Hits containers
    art::Handle< std::vector<crt::CRTHit> > crtHitHandle;  // Container of crt hits
    std::vector<art::Ptr<crt::CRTHit> >     crtHits;       // Vector of crt hits 
    // if there's none, kill the event
    if(!e.getByLabel(fCRTHitLabel, crtHitHandle)) {  e.put(std::move(VtxCol)); e.put(std::move(vertex_track_assn_v)); return false;}
    // fill the container of crt hits                                                                                                                        
    art::fill_ptr_vector(crtHits, crtHitHandle);
    // Find associations
    art::FindOneP<recob::Track> fCRT2TPC(crtHitHandle, e , fCRTTrackAssnProducer);
    if (fCRT2TPC.isValid())
      {
	for (unsigned int indexAssn = 0; indexAssn < fCRT2TPC.size(); ++indexAssn )
	  {
	    cet::maybe_ref<recob::Track const> trackCRT2TPC(*fCRT2TPC.at(indexAssn));
	    if (!trackCRT2TPC) continue;
	    recob::Track const& aTrack(trackCRT2TPC.ref());
	    matchedRecoTrkKey.push_back(aTrack.ID()); 
	  }// End Loop on Association
      }// if there's an associated track                                                                                                                        
    // Now we know what ID identifies the reco track we're interested in.                                                                                       
  }// if use CRT information
  //if (recoTrk->ID() != matchedRecoTrkKey) continue;
  




  // if this track associated with a crt hit?

  /*Scope of this bit of code is to create an association between bad vtx and tracks */
  // Let's determe if our verteces are good or not
  bool thisVtxIsGood = true;
  //Loop on verteces

  _nVtxs         = vtxlist.size();
  _nTracks       = tracklist.size();
  _nTrackCRTAssn = matchedRecoTrkKey.size();
  for (size_t iVtx = 0; iVtx < vtxlist.size(); iVtx++)
    {
      // shortest distance between this vtx and the considered tracks
      float shortestDistance = 99999.;
      thisVtxIsGood = true;
      auto thisVtx = vtxlist[iVtx];
      auto vPos = thisVtx->position();
      if (fVerbose) std::cout<<"__ Vtx Pos: "<< vPos.X()<<" "<<vPos.Y()<<" "<<vPos.Z() <<"__\n";
      int trackUsed = 0;
      for ( auto const& thisTrack : tracklist )	{
	if (fUseCRTTag) { // this is where we check that the track has a crt assn
	  int trackID = thisTrack->ID();
	  // I need to check if this track ID is in the keys
	  // if it is, I'm interested in this track. Otherwise, I'm not and I shall continue.
	  if(! (std::find(matchedRecoTrkKey.begin(), matchedRecoTrkKey.end(), trackID) != matchedRecoTrkKey.end())) continue;
	}// are we using the CRT info?
	trackUsed++;
	// Grab this track TrajectorPoints & calculate the trajectory
	for(size_t p = 1; p < thisTrack->NumberTrajectoryPoints(); ++p){
	  const auto& pos_cur = thisTrack->LocationAtPoint(p);	  
	  float dx = vPos.X() - pos_cur.x(); 
	  float dy = vPos.Y() - pos_cur.y(); 
	  float dz = vPos.Z() - pos_cur.z();
	  
	  float thisTrajPtDistance = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
	  if (fVerbose){
	    std::cout<<" dx "<< vPos.X() <<" "<<pos_cur.x() <<"\n"; 
	    std::cout<<" dy "<< vPos.Y() <<" "<<pos_cur.y() <<"\n"; 
	    std::cout<<" dz "<< vPos.Z() <<" "<<pos_cur.z() <<"\n";
	    std::cout<<" distance "<< thisTrajPtDistance <<"\n";
	  }
	  // calculate the shortest distance for this vtx
	  if (shortestDistance > thisTrajPtDistance ) shortestDistance = thisTrajPtDistance;
	  // if the distance between the space point and the vertex is less than the minimum, 
	  // the vtx is not good. Stop looping on the space points (and on the tracks later on)
	  if (thisTrajPtDistance < fMinDistance) {thisVtxIsGood = false; break;} 
	}
	_closestDistanceVtxTrk.push_back(shortestDistance);
	// If the vtx is already busted, save track-vtx association
	// stop looping on the tracks, move on to the new vtx!
	if (thisVtxIsGood == false) { vertex_track_assn_v->addSingle(thisVtx, thisTrack); break;} 
      } // The loop on the tracks
      // if after I looped on all the tracks the vertex is still good, save the vertex in the event
      if (thisVtxIsGood == true) {  VtxCol->push_back(*thisVtx);} 
    }// The loop on the vtx


  _nGoodVtx = (*VtxCol).size() ;
  if ((*VtxCol).size() > 0 ) _keepEvent = true;
  // Fill the tree.
  _tree->Fill();

  // Add the data products to the event.
  e.put(std::move(VtxCol));
  e.put(std::move(vertex_track_assn_v));
  
  // Return false if you want to discart the event
  if (fuseAsFilter && _keepEvent == false)
    return false;

  // Otherwise return true
  return true;

} // End of the filter module


DEFINE_ART_MODULE(CRTDistanceFilter)

//  LocalWords:  XYZ
