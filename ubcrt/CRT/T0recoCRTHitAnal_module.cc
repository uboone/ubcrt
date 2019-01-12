////////////////////////////////////////////////////////////////////////
// Class:       T0recoCRTHitAnal
// Module Type: analyzer
// File:        T0recoCRTHitAnal_module.cc
//
//   Written Dec 2018 to test developments in matching TPC tracks
//      to CRT hits
//
//   Last updated 12/19/2018
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

//data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/Utilities/AssociationUtil.h"

//CRT data-products                                                              
#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTTrack.hh"
#include "ubcrt/CRT/CRTAuxFunctions.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

//Root
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3S.h"
#include "TProfile.h"
#include "TF1.h"
#include "TDatime.h"
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <typeinfo>

#include <memory>

namespace crt {
  class T0recoCRTHitAnal;
}

class crt::T0recoCRTHitAnal : public art::EDAnalyzer {
public:
  explicit T0recoCRTHitAnal(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  T0recoCRTHitAnal(T0recoCRTHitAnal const &) = delete;
  T0recoCRTHitAnal(T0recoCRTHitAnal &&) = delete;
  T0recoCRTHitAnal & operator = (T0recoCRTHitAnal const &) = delete;
  T0recoCRTHitAnal & operator = (T0recoCRTHitAnal &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void SortTrackPoints (const recob::Track& track, std::vector<TVector3>& sorted_trk);

private:

  // Declare member data here.
  art::ServiceHandle<art::TFileService> tfs;

  std::string  data_label_TPCTrack;
  std::string  data_label_crtT0;
  std::string  data_label_acptT0;
  float fTQCutOpAng;
  float fTQCutLength;
  bool fverbose;
  bool fIsMC;
  
  TH1F* hDiffT_CRT_Flash;
  TH1F* hDiffT_CRT_AFlash;
  TH1F* hGeoMatch;
  TH1F* hTheta;
  TH1F* hPhi;
  TH1F* hLength;
  TH1F* hOpAng;
  TH1F* hMTheta;
  TH1F* hMPhi;
  TH1F* hMTime;
  TH1F* hMLength;
  TH1F* hMOpAng;
  TH1F* hMlowx;
  TH1F* hMhighx;
  TH1F* hATheta;
  TH1F* hAPhi;
  TH1F* hATime;
  TH1F* hALength;
  TH1F* hAOpAng;
  TH1F* hAlowx;
  TH1F* hAhighx;
  TH1F* hAMTheta;
  TH1F* hAMPhi;
  TH1F* hAMTime;
  TH1F* hAMLength;
  TH1F* hAMOpAng;
  TH1F* hAMlowx;
  TH1F* hAMhighx;

};


crt::T0recoCRTHitAnal::T0recoCRTHitAnal(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  // data_labelCRThit_(p.get<std::string>("data_labelCRThit")),
  // data_label_flash_(p.get<std::string>("data_label_flash_")),
  // data_label_DAQHeader_(p.get<std::string>("data_label_DAQHeader_")),
  data_label_TPCTrack(p.get<std::string>("data_label_TPCtrack","pandoraCosmic")),
  data_label_crtT0(p.get<std::string>("data_label_crtT0","t0recocrthit")),
  data_label_acptT0(p.get<std::string>("data_label_acptT0","pandoraComsicT0Reco")),
  fTQCutOpAng(p.get<float>("TQCutOpAng",0.95)),
  fTQCutLength(p.get<float>("TQCutLength",20)),
  fverbose(p.get<bool>("verbose",false)),
  fIsMC(p.get<bool>("IsMC",false))
 // More initializers here.
{}

void crt::T0recoCRTHitAnal::analyze(art::Event const & evt)
{


  if (fIsMC) std::cout << "this is monte carlo" << std::endl;
  
  //get TPC Tracks  
  art::Handle< std::vector<recob::Track> > rawHandle_TPCtrack;
  evt.getByLabel(data_label_TPCTrack, rawHandle_TPCtrack);
  //check to make sure the data we asked for is valid                                                                                                         
  if(!rawHandle_TPCtrack.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " recob::Track " << " in module " << data_label_TPCTrack << std::endl;
    std::cout << std::endl;
    return;
  }
  //get better access to the data	    
  
  std::vector<recob::Track> const& TPCTrackCollection(*rawHandle_TPCtrack);
  if(fverbose){
    std::cout<<"  TPCTrackCollection.size()  "<<TPCTrackCollection.size()<<std::endl;
    
  }
  //get TPCTracks                                                                                                                                               

  //get Optical Flash			   
  
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlash;
  evt.getByLabel(data_label_flash_, rawHandle_OpFlash);
  std::vector<recob::OpFlash> const& OpFlashCollection(*rawHandle_OpFlash);
  if(fverbose){
    std::cout<<"  OpFlashCollection.size()  "<<OpFlashCollection.size()<<std::endl;
  }
  //get Optical Flash 
  
                 
  // look for match between CRT and flashes   
  for(std::vector<int>::size_type i = 0; i != OpFlashCollection.size(); i++) {//B
      
    recob::OpFlash my_OpFlash = OpFlashCollection[i];
      
    auto Yflash = my_OpFlash.YCenter();
    auto Zflash = my_OpFlash.ZCenter();
    auto PEflash = my_OpFlash.TotalPE();
    auto Timeflash = my_OpFlash.Time(); //in us from trigger time
    auto Timeflash_ns = (Timeflash * 1000);
    auto Timeflash_ns_GPS = evt_timeGPS_nsec + (Timeflash * 1000);      
    int fbeam = my_OpFlash.OnBeamTime();
    uint32_t Flash_sec = evt_timeGPS_sec;
      
    //    hFlashTimeDis->Fill(Timeflash);
    
    //loop over CRT tzeros
    float min_deltat = 3000.0; int best_time_match = -1;
    if (tzerolist.size()>0) {   
      for(size_t tzIter = 0; tzIter < tzerolist.size(); ++tzIter){   
	float diff;
	if (fTimeSelect==0) diff = fabs(0.001*(tzerolist[tzIter]->ts0_ns+fTimeZeroOffset-(double)evt_timeGPS_nsec)-Timeflash);
	else diff= fabs(0.001*(tzerolist[tzIter]->ts1_ns+fHardDelay)-Timeflash);
	if (diff<min_deltat) { min_deltat=diff; best_time_match=tzIter;}
      } // loop over tzeros     
      if (best_time_match>=0) {
	hDiffT_CRT_Flash->Fill(min_deltat);	
	if (fverbose) std::cout << "Closest CRT hit in time to this flash is tzero no " << 
			best_time_match << " at time (us) " <<
			0.001*(tzerolist[best_time_match]->ts1_ns+fHardDelay) << std::endl;
      }// if CRT-flash match found
    } // if tzeros 

    if(fverbose){ 
      std::cout<<"event: "<<fEvtNum<<std::endl;
      std::cout<<"Flash: "<<i<<std::endl;
      std::cout<<"Beam: "<<fbeam<<std::endl;
      std::cout<<"Zflash: "<<Zflash<<std::endl;
      std::cout<<"Yflash: "<<Yflash<<std::endl;
      std::cout<<"PEflash: "<<PEflash<<std::endl;
      std::cout.precision(19);
      std::cout<<" "<<std::endl;
      std::cout<<"Flash time: "<<Flash_sec<< " seconds"<<std::endl;
      std::cout<<"Flash time: "<<Timeflash<< "  us w.r.t to trigger"<<std::endl;
      std::cout<<"Flash time: "<<Timeflash_ns<< "  ns w.r.t to trigger"<<std::endl;
      std::cout<<"Flash time: "<<Timeflash_ns_GPS<< "  ns in GPS"<<std::endl;
      std::cout<<" "<<std::endl;
    }
  } // loop over flashes
   

  // 1.11436 mm/us   
  double driftvel = 0.11436; //   units  cm/us 
  // const detinfo::DetectorProperties *_detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // double const vdrift =  _detprop->DriftVelocity();  

  // grab T0 objects associated with tracks    
  art::FindMany<anab::T0> trk_t0_assn_v(rawHandle_TPCtrack, evt, data_label_crtT0);
  // grab flashes associated with tracks (anode or cathode crossers)
  art::FindMany<anab::T0> trk_flash_assn_v(rawHandle_TPCtrack, evt, data_label_acptT0 );
  
  for(std::vector<int>::size_type i = 0; i != TPCTrackCollection.size(); i++) {     
    recob::Track my_TPCTrack = TPCTrackCollection[i];    
    const std::vector<const anab::T0*>& T0_v = trk_t0_assn_v.at(i);
    const std::vector<const anab::T0*>& T0_acpt = trk_flash_assn_v.at(i);
    
    double t_theta=my_TPCTrack.Theta();
    double t_phi=my_TPCTrack.Phi();
    if (t_phi>0) {t_phi-=3.14159; t_theta=3.14159-my_TPCTrack.Theta();}
    double t_len = my_TPCTrack.Length();
    // get track directional cosines
    double trackCosStart[3]={0.,0.,0.};
    double trackCosEnd[3]={0.,0.,0.};
    my_TPCTrack.Direction(trackCosStart,trackCosEnd);      
    double t_opang = trackCosStart[0]*trackCosEnd[0] +  trackCosStart[1]*trackCosEnd[1] + 
      trackCosStart[2]*trackCosEnd[2];
    if (t_len>fTQCutLength && t_opang>fTQCutOpAng) {
      //    if (t_len>20 && t_opang>0.95) {
      auto startP = my_TPCTrack.Start();
      auto endP = my_TPCTrack.End();
      double lowx = startP.X();
      double highx = endP.X();
      if (lowx>highx) {
	highx=lowx;
	lowx=endP.X();
      }
      
      hTheta->Fill(t_theta);
      hPhi->Fill(t_phi);
      hOpAng->Fill(t_opang);
      hLength->Fill(t_len);
      
      double atracktime=0.0;
      bool b_acpt = false;
      if (T0_acpt.size()==1) { 
	b_acpt=true;
	auto t0 = T0_acpt.at(0);
	atracktime = t0->Time();  //track time in us, t0->time() in us
	if (fverbose) std::cout << "acpt T0 time and trigger type  " << 
		      atracktime << " " << t0->TriggerType() << std::endl;      
	hATheta->Fill(t_theta);
	hAPhi->Fill(t_phi);
	hATime->Fill(atracktime);
	hAOpAng->Fill(t_opang);
	hALength->Fill(t_len);
	hAlowx->Fill(lowx-atracktime*driftvel);
	hAhighx->Fill(highx-atracktime*driftvel);
      }

      double ctracktime=0.0;
      if (T0_v.size() == 1){
	auto t0 = T0_v.at(0);
	ctracktime = t0->Time();  //track time in us, t0->time() in us
	if (fverbose) std::cout << "CRT T0 time and trigger type  " << 
			ctracktime << " " << t0->TriggerType() << std::endl;
	hMTheta->Fill(t_theta);
	hMPhi->Fill(t_phi);
	hMTime->Fill(ctracktime);
	hMOpAng->Fill(t_opang);
	hMlowx->Fill(lowx-ctracktime*driftvel);
	hMhighx->Fill(highx-ctracktime*driftvel);
	hMLength->Fill(t_len);
	if (b_acpt) {
	  hAMTheta->Fill(t_theta);
	  hAMPhi->Fill(t_phi);
	  hDiffT_CRT_AFlash->Fill(atracktime-ctracktime);
	  hAMTime->Fill(ctracktime);
	  hAMOpAng->Fill(t_opang);
	  hAMLength->Fill(t_len);
	  hAMlowx->Fill(lowx-atracktime*driftvel);
	  hAMhighx->Fill(highx-atracktime*driftvel);
	}
      }
    }// if (t_len>5 && t_opang>0.8)
      
  }//loop over tracks
    
    
    
}

void crt::T0recoCRTHitAnal::beginJob()
{

  hDiffT_CRT_AFlash = tfs->make<TH1F>("hDiffT_CRT_AFlash","hDiffT_CRT_AFlash",400,-2,2);
  hDiffT_CRT_AFlash->GetXaxis()->SetTitle("Flash_time_ACPT - CRTHit_Time (us)");
  hDiffT_CRT_AFlash->GetYaxis()->SetTitle("Entries/bin");

  hDiffT_CRT_Flash = tfs->make<TH1F>("hDiffT_CRT_Flash","hDiffT_CRT_Flash",400,-2,2);
  hDiffT_CRT_Flash->GetXaxis()->SetTitle("Nearest_Flash - CRTHit_Time (us)");
  hDiffT_CRT_Flash->GetYaxis()->SetTitle("Entries/bin");

  // hDiffT_T0Flash = tfs->make<TH1F>("hDiffT_T0Flash","hDiffT_T0Flash",5000,-20000,20000);
  // hDiffT_T0Flash->GetXaxis()->SetTitle("T0Time - T0Time_Flash (ns)");
  // hDiffT_T0Flash->GetYaxis()->SetTitle("Entries/bin");

  // hGeoMatch = tfs->make<TH1F>("hGeoMatch","hGeoMatch",5000,-20000,20000);
  // hGeoMatch->GetXaxis()->SetTitle("FlashTime - CRTTrack_Time (ns)");
  // hGeoMatch->GetYaxis()->SetTitle("Entries/bin");

  hTheta = tfs->make<TH1F>("hTheta","hTheta",30,0,3.14159);
  hPhi = tfs->make<TH1F>("hPhi","hPhi",30,-3.14159,0);
  hMTheta = tfs->make<TH1F>("hMTheta","hMTheta",30,0,3.14159);
  hMPhi = tfs->make<TH1F>("hMPhi","hMPhi",30,-3.14159,0);
  hATheta = tfs->make<TH1F>("hATheta","hATheta",30,0,3.14159);
  hAPhi = tfs->make<TH1F>("hAPhi","hAPhi",30,-3.14159,0);
  hAMTheta = tfs->make<TH1F>("hAMTheta","hAMTheta",30,0,3.14159);
  hAMPhi = tfs->make<TH1F>("hAMPhi","hAMPhi",30,-3.14159,0);
  hMTime = tfs->make<TH1F>("hMtime","htime",200,-4000,4000);
  hATime = tfs->make<TH1F>("hAMtime","hAtime",200,-4000,4000);
  hAMTime = tfs->make<TH1F>("hAMtime","hAtime",200,-4000,4000);
  hLength = tfs->make<TH1F>("hLength","hLength",80,0,400.0);
  hLength->GetXaxis()->SetTitle("Track Length (cm)");
  hOpAng = tfs->make<TH1F>("hOpAng","hOpAng",60,0.85,1.00);
  hOpAng->GetXaxis()->SetTitle("cos(alpha)");
  hMLength = tfs->make<TH1F>("hMLength","hMLength",80,0,400.0);
  hMLength->GetXaxis()->SetTitle("Track Length (cm)");
  hMOpAng = tfs->make<TH1F>("hMOpAng","hMOpAng",60,0.85,1.00);
  hMOpAng->GetXaxis()->SetTitle("cos(alpha)");
  hMlowx = tfs->make<TH1F>("hMlowx","hMlowx",300,-10.,290.);
  hMlowx->GetXaxis()->SetTitle("x (cm)");
  hMhighx = tfs->make<TH1F>("hMhighx","hMhighx",300,-10.,290.);
  hMhighx->GetXaxis()->SetTitle("x (cm)");
  hALength = tfs->make<TH1F>("hALength","hALength",80,0,400.0);
  hALength->GetXaxis()->SetTitle("Track Length (cm)");
  hAOpAng = tfs->make<TH1F>("hAOpAng","hAOpAng",60,0.85,1.00);
  hAOpAng->GetXaxis()->SetTitle("cos(alpha)");
  hAlowx = tfs->make<TH1F>("hAlowx","hAlowx",300,-10.,290.);
  hAlowx->GetXaxis()->SetTitle("x (cm)");
  hAhighx = tfs->make<TH1F>("hAhighx","hAhighx",300,-10.,290.);
  hAhighx->GetXaxis()->SetTitle("x (cm)");
  hAMLength = tfs->make<TH1F>("hAMLength","hAMLength",80,0,400.0);
  hAMLength->GetXaxis()->SetTitle("Track Length (cm)");
  hAMOpAng = tfs->make<TH1F>("hAMOpAng","hAMOpAng",60,0.85,1.00);
  hAMOpAng->GetXaxis()->SetTitle("cos(alpha)");
  hAMlowx = tfs->make<TH1F>("hAMlowx","hAMlowx",300,-10.,290.);
  hAMlowx->GetXaxis()->SetTitle("x (cm)");
  hAMhighx = tfs->make<TH1F>("hAMhighx","hAMhighx",300,-10.,290.);
  hAMhighx->GetXaxis()->SetTitle("x (cm)");

}

void crt::T0recoCRTHitAnal::endJob()
{
  // Implementation of optional member function here.

  // Double_t norm = hFlashCountClone->GetEntries();
  // hFlashCountClone->Scale(1/norm);


}


DEFINE_ART_MODULE(crt::T0recoCRTHitAnal)