#ifndef crt_CRTAnaFun_hh
#define crt_CRTAnaFun_hh


#include <stdio.h>
#include <map> 
#include <vector> 

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include "cetlib/search_path.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

//#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//#include "art/Framework/Services/Optional/TFileService.h"
#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTTrack.hh"
#include "ubcrt/CRT/CRTAuxFunctions.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"


#include "lardataobj/RecoBase/Track.h"                                                                
#include "lardataobj/RecoBase/Hit.h"                                                                  
#include "lardataobj/AnalysisBase/T0.h"                                                               
#include "lardataobj/AnalysisBase/CosmicTag.h"                                                        
#include "lardataobj/AnalysisBase/Calorimetry.h"                                                      
#include "lardataobj/MCBase/MCTrack.h"                                                                
#include "lardataobj/RecoBase/OpFlash.h"                                                              
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

//#include "Event_Tree.h"


#include "Pandora/PdgTable.h"


namespace crtana {

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;
    typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
    typedef std::vector< art::Ptr<recob::Track> > TrackVector;
    typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;


    

  namespace auxfunc {
  
    

    int test(int a);
    //int test2(int a);
    //std::vector<crt::CRTHit> getCRTHitVector(art::Event const & evt, std::string data_label_hits_);
    void GetPFParticleIdMap(const crtana::PFParticleHandle &pfParticleHandle, crtana::PFParticleIdMap &pfParticleMap);
    void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles);
    void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,std::string m_trackLabel, std::string m_showerLabel);
    
    double TpcTrack_match_CrtHit(const art::Ptr<recob::Track> &ptracks, std::vector<crt::CRTHit> const& CRTHitCollection, std::vector<double> &dist, double crthitmatch_);
    void SortTrackPoints(const recob::Track& track, std::vector<TVector3>& sorted_trk);


  }

}

#endif
