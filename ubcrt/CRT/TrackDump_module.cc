////////////////////////////////////////////////////////////////////////
// Class:       TrackDump
// Module Type: analyzer
// File:        TrackDump_module.cc
//
// Generated at Mon Jul  3 03:51:03 2017 by David Lorca Galindo using artmod
// from cetpkgsupport v1_11_00.
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

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <artdaq-core/Data/Fragment.hh>

#include "art/Framework/Services/Optional/TFileService.h"

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTTzero.hh"
#include "ubobj/CRT/CRTTrack.hh"
#include "ubcrt/CRT/CRTAuxFunctions.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

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


const int kMaxCRThits = 1000;
const int kMaxCRTtzs = 1000;
const int kMaxCRTtracks = 1000;
const int kMaxTPCtracks = 100;
const int kMaxPMTflash = 100;


 // namespace crt {
 //   class TrackDump;
 // }

class TrackDump : public art::EDAnalyzer {
public:
  explicit TrackDump(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackDump(TrackDump const &) = delete;
  TrackDump(TrackDump &&) = delete;
  TrackDump & operator = (TrackDump const &) = delete;
  TrackDump & operator = (TrackDump &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  void ResetVars();

  art::ServiceHandle<art::TFileService> tfs;
  // Declare member data here.
  
  uint32_t fEvtNum; //Number of current event                       
  uint32_t frunNum;                //Run Number taken from event  
  uint32_t fsubRunNum;             //Subrun Number taken from event         
  std::string  fTrackModuleLabel;
  bool fSaveTPCTrackInfo;
  bool fSavePMTFlashInfo;
  std::string  data_labeltrack_;
  std::string  data_labeltzero_;
  std::string  data_label_t0A_;
  std::string  data_label_t0C_;
  std::string  data_labelhit_;
  std::string  data_label_flash_;
  std::string  data_label_DAQHeader_;
  //  bool fIsMC;
  int fHardDelay_;
  int fTimeZeroOffset;
  int fDriftVel;
  int verbose_;
  


  //art::InputTag opFlashTag("opflashSat");


  //quality plots

  TH2F* hplavspla;
  TH1F* hTlength;
  TH1F* hTtime;
  TH2F* hTlengthvsTime;
  TH2F* hTlengthvsTimeAbs;
  TProfile* hTlengthvsTimeAbs_prof;
  TH1F* htheta;
  TH1F* hphi;
  TH1F* hts0_ns;
  TH2F* hTvsH;

  TH2F* HitDistBot;
  TH2F* HitDistFT;
  TH2F* HitDistPipe;
  TH2F* HitDistTop;
 //quality plots                                                                                                                                          

  TTree*       fTree;
  // run information
  int run;
  int subrun;
  int event;
  double evttime;
  // CRT hits
  int nCRThits;
  int hit_plane[kMaxCRThits];
  double hit_time_s[kMaxCRThits];
  double hit_time0[kMaxCRThits];
  double hit_time1[kMaxCRThits];
  double hit_charge[kMaxCRThits];
  double hit_posx[kMaxCRThits];
  double hit_posy[kMaxCRThits];
  double hit_posz[kMaxCRThits]; 
  int hit_feb1[kMaxCRThits]; 
  int hit_strip1[kMaxCRThits]; 
  int hit_feb2[kMaxCRThits]; 
  int hit_strip2[kMaxCRThits]; 
  double hit_pe1[kMaxCRThits]; 
  double hit_pe2[kMaxCRThits]; 
  double hit_sipm1a[kMaxCRThits]; 
  double hit_sipm1b[kMaxCRThits]; 
  double hit_sipm2a[kMaxCRThits]; 
  double hit_sipm2b[kMaxCRThits]; 
  int hit_mcflag[kMaxCRThits]; 
  // CRT Tzeros
  int nCRTtzeros;
  double tz_time_s[kMaxCRTtzs];
  double tz_time0[kMaxCRTtzs];
  double tz_time1[kMaxCRTtzs];
  int tz_hits0[kMaxCRTtzs];
  int tz_hits1[kMaxCRTtzs];
  int tz_hits2[kMaxCRTtzs];
  int tz_hits3[kMaxCRTtzs];
  double tz_pes0[kMaxCRTtzs];
  double tz_pes1[kMaxCRTtzs];
  double tz_pes2[kMaxCRTtzs];
  double tz_pes3[kMaxCRTtzs];
  // CRT tracks
  int nCRTtracks;
  double ct_theta[kMaxCRTtracks];
  double ct_phi[kMaxCRTtracks];
  double ct_length[kMaxCRTtracks];
  double ct_time_sec[kMaxCRTtracks];
  double ct_time0[kMaxCRTtracks];
  double ct_time1[kMaxCRTtracks];
  double ct_x1[kMaxCRTtracks];
  double ct_y1[kMaxCRTtracks];
  double ct_z1[kMaxCRTtracks];
  double ct_x2[kMaxCRTtracks];
  double ct_y2[kMaxCRTtracks];
  double ct_z2[kMaxCRTtracks];
  // TPC tracks
  int nTPCtracks;
  double trkstartx[kMaxTPCtracks];
  double trkstarty[kMaxTPCtracks];
  double trkstartz[kMaxTPCtracks];
  double trkendx[kMaxTPCtracks];
  double trkendy[kMaxTPCtracks];
  double trkendz[kMaxTPCtracks];
  double trkstartdcosx[kMaxTPCtracks];
  double trkstartdcosy[kMaxTPCtracks];
  double trkstartdcosz[kMaxTPCtracks];
  double trkenddcosx[kMaxTPCtracks];
  double trkenddcosy[kMaxTPCtracks];
  double trkenddcosz[kMaxTPCtracks];
  double trktheta[kMaxTPCtracks];
  double trkphi[kMaxTPCtracks];
  double trklen[kMaxTPCtracks];
  double tzeroACPT[kMaxTPCtracks];
  double tzeroCRT[kMaxTPCtracks];
  //  double dcaCRT[kMaxTPCtracks];
  int planeCRT[kMaxTPCtracks];
  int feb1CRT[kMaxTPCtracks];
  int feb2CRT[kMaxTPCtracks];
  double pes1CRT[kMaxTPCtracks];
  double pes2CRT[kMaxTPCtracks];
  double xCRT[kMaxTPCtracks];
  double yCRT[kMaxTPCtracks];
  double zCRT[kMaxTPCtracks];
  double angCRT[kMaxTPCtracks];
  // PMT flash
  int nPMTflash;
  double fl_time[kMaxPMTflash];
  double fl_pe[kMaxPMTflash];
  double fl_y[kMaxPMTflash];
  double fl_z[kMaxPMTflash];

  
};


TrackDump::TrackDump(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    fTrackModuleLabel(p.get<std::string>("TrackModuleLabel")),
    fSaveTPCTrackInfo(p.get< bool >("SaveTPCTrackInfo", false)), 
    fSavePMTFlashInfo(p.get< bool >("SavePMTFlashInfo", false)), 
    data_labeltrack_(p.get<std::string>("data_labeltrack")),
    data_labeltzero_(p.get<std::string>("data_labeltzero")),
    data_label_t0A_(p.get<std::string>("data_label_t0ACPT", "pandoraCosmicT0Reco" )),
    data_label_t0C_(p.get<std::string>("data_label_t0CRT")),
    data_labelhit_(p.get<std::string>("data_labelhit")),
    data_label_flash_(p.get<std::string>("data_label_flash_")),
    data_label_DAQHeader_(p.get<std::string>("data_label_DAQHeader_")),
    fHardDelay_(p.get<int>("HardDelay",40000)),
    fTimeZeroOffset(p.get<int>("TimeZeroOffset",60000)),
    fDriftVel(p.get<float>("DriftVel",0.111436)),   // cm/us
    verbose_(p.get<int>("verbose"))
    // More initializers here.    
{
}

void TrackDump::analyze(art::Event const & evt)
{
  ResetVars();

  frunNum    = evt.run();
  fsubRunNum = evt.subRun();
  fEvtNum = evt.event();
  
  art::Timestamp evtTime = evt.time();
  auto evt_time_sec = evtTime.timeHigh();
  auto evt_time_nsec = evtTime.timeLow();
  double evt_timeGPS_sec = 0.0;
  double evt_timeGPS_nsec = 0.0;
  double evt_timeNTP_sec = 0.0;
  double evt_timeNTP_nsec = 0.0;
  double timstp_diff = 0.0;

  //  fIsMC=1;
  if (evt.isRealData()) {  // get DAQ timestamps if this is data
    //    fIsMC=0;
   //get DAQ Header                                                                  
  //Commentar para old swizzler, sin DAQ Header
    art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;  
    evt.getByLabel(data_label_DAQHeader_, rawHandle_DAQHeader);
    
    //check to make sure the data we asked for is valid 
    if(!rawHandle_DAQHeader.isValid()){
      std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
		<< ", event " << evt.event() << " has zero"
		<< " DAQHeaderTimeUBooNE  " << " in with label " << data_label_DAQHeader_ << std::endl;    
      return;
    }
    
    
    raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
    art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();  
    evt_timeGPS_sec = evtTimeGPS.timeHigh();
    evt_timeGPS_nsec = (double)evtTimeGPS.timeLow();
    art::Timestamp evtTimeNTP = my_DAQHeader.ntp_time();
    evt_timeNTP_sec = evtTimeNTP.timeHigh();
    evt_timeNTP_nsec = (double)evtTimeNTP.timeLow();
    timstp_diff = std::abs(evt_timeGPS_nsec - evt_timeNTP_nsec);
    
    if(verbose_==1){
      std::cout<< "Run:  "<<frunNum << "   subRun: " <<fsubRunNum<<std::endl;
      std::cout<<"event: "<<fEvtNum <<std::endl;
      std::cout.precision(19);
      std::cout<<"  GPS time second:  "<<evt_timeGPS_sec<<std::endl;
      std::cout<<"  GPS time nano_second:  "<<evt_timeGPS_nsec<<std::endl;
      std::cout<<"  NTP time second:  "<<evt_timeNTP_sec<<std::endl;    
      std::cout<<"  NTP time nano_second:  "<<evt_timeNTP_nsec<<std::endl;
      std::cout<<"  event time second:  "<<evt_time_sec<<std::endl;
      std::cout<<"  event time nano_second:  "<<evt_time_nsec<<std::endl;
      std::cout<<"  difference between GPS and NTP:  "<<evt_timeGPS_nsec - evt_timeNTP_nsec<<" ns"<<std::endl;
      std::cout<<"  ABS difference between GPS and NTP:  "<<timstp_diff<<" ns"<<std::endl;
      
      if( (evt_time_sec==evt_timeGPS_sec) && (evt_time_nsec==evt_timeGPS_nsec))  std::cout<<" Event time type is: GPS  "<<std::endl;
      if( (evt_time_sec==evt_timeNTP_sec) && (evt_time_nsec==evt_timeNTP_nsec))  std::cout<<" Event time type is: NTP  "<<std::endl;
      //getchar();
    }  
  }
  

  if (fSaveTPCTrackInfo) {
    
      // get TPC Track List 
    art::Handle< std::vector<recob::Track>  > trackListHandle; 
    std::vector<art::Ptr<recob::Track> >  tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
      art::fill_ptr_vector(tracklist, trackListHandle);
    //check to make sure the data we asked for is valid
    if(!trackListHandle.isValid()){
      std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
		<< ", event " << evt.event() << " has zero"
		<< " tracks " << " in module " << fTrackModuleLabel << std::endl;
      std::cout << std::endl;
      return;
    }
  //check whether tzeros exist
    bool iT0acpt=false;
    bool iT0crt = false;
    art::Handle< std::vector<anab::T0> > rawHandle_Tzero;
    evt.getByLabel(data_label_t0A_, rawHandle_Tzero);
    if(rawHandle_Tzero.isValid()) { iT0acpt=true;
      // grab flashes associated with tracks (anode or cathode crossers)
      //      art::FindMany<anab::T0> trk_t0A_assn_v(trackListHandle, evt, data_label_t0A_);
      //      std::cout << "found data product for acpt times" << std::endl;
    }
    evt.getByLabel(data_label_t0C_, rawHandle_Tzero);
    if(rawHandle_Tzero.isValid()) {iT0crt=true;
      // grab T0 objects associated with tracks    
      //      art::FindMany<anab::T0> trk_t0C_assn_v(trackListHandle, evt, data_label_t0C_);
      std::cout << "found data product for crt times" << std::endl;
    }
    auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
    
    double const vdrift = fDriftVel;  // in cm/us
    
    nTPCtracks = tracklist.size();
    if (nTPCtracks>kMaxTPCtracks) nTPCtracks=kMaxTPCtracks;
    for(int j = 0; j < nTPCtracks; j++) {
      if (verbose_)std::cout<<"in iTrackloop\n";
      art::Ptr<recob::Track> ptrack(trackListHandle, j);
      const recob::Track& track = *ptrack;
      
      auto pos       = track.Vertex();
      auto dir_start = track.VertexDirection();
      auto dir_end   = track.EndDirection();
      auto startP       = track.Start();
      auto endP       = track.End();

      //
      trklen[j]= track.Length(); //length(track);
      trkstartx[j]=pos.X();
      trkstarty[j]=pos.Y();
      trkstartz[j]=pos.Z();
      trkendx[j]=endP.X();
      trkendy[j]=endP.Y();
      trkendz[j]=endP.Z();
      trkstartdcosx[j]=dir_start.X();
      trkstartdcosy[j]=dir_start.Y();
      trkstartdcosz[j]=dir_start.Z();
      trkenddcosx[j]=dir_end.X();
      trkenddcosy[j]=dir_end.Y();
      trkenddcosz[j]=dir_end.Z();
      trktheta[j]=dir_start.Theta();
      trkphi[j]=dir_start.Phi();
      tzeroACPT[j]=-9999.0;
      if (iT0acpt) { 
	art::FindMany<anab::T0> trk_t0A_assn_v(trackListHandle, evt, data_label_t0A_);
	const std::vector<const anab::T0*>& T0_acpt = trk_t0A_assn_v.at(j);
	if (T0_acpt.size()==1) { 
	  auto t0 = T0_acpt.at(0);
	  tzeroACPT[j]=t0->Time();  //track time in us, t0->time() in us
	}
      }
      tzeroCRT[j]=-9999.0; planeCRT[j]=-1;
      if (iT0crt) {
	if (verbose_)std::cout<<"in iT0crt\n";
	//	art::FindMany<anab::T0> trk_t0C_assn_v("pandora", evt, "pandoraCrtHitMatch");
	art::FindMany<anab::T0> trk_t0C_assn_v(trackListHandle, evt, data_label_t0C_);
	const std::vector<const anab::T0*>& T0_v = trk_t0C_assn_v.at(j);
	if (T0_v.size()==1) { 
	  auto t0 = T0_v.at(0);
	  tzeroCRT[j]=t0->Time();	  
	  planeCRT[j]=t0->TriggerBits();  
	  //	  std::cout << "dca value is " << t0->TriggerConfidence() << std::endl;
	  art::FindMany<crt::CRTHit> trk_hit_assn_v(trackListHandle, evt, data_label_t0C_);
	  const std::vector<const crt::CRTHit*>& hit_v = trk_hit_assn_v.at(j);
	  if (hit_v.size()>0) {
	    auto thishit = hit_v.at(0);
	    feb1CRT[j]=thishit->feb_id[0];	    feb2CRT[j]=thishit->feb_id[1];
	    xCRT[j]=thishit->x_pos;	    yCRT[j]=thishit->y_pos;	    zCRT[j]=thishit->z_pos;

	    std::vector<std::pair<int,float>> pes1 = thishit->pesmap.find(feb1CRT[j])->second; 
	    std::vector<std::pair<int,float>> pes2 = thishit->pesmap.find(feb2CRT[j])->second; 
	    if (pes1.size()==2) { // is simulated hit
	      std::pair<int,float> ind_pes1 = pes1[0];
	      std::pair<int,float> ind_pes2 = pes1[1];
	      pes1CRT[j]=ind_pes1.second+ind_pes2.second;
	      std::pair<int,float> ind2_pes1 = pes2[0];
	      std::pair<int,float> ind2_pes2 = pes2[1];
	      pes2CRT[j]=ind2_pes1.second+ind2_pes2.second;
	    }
	    else {
	      float hitpe1 = 0; 
	      for (uint isp=0;isp<pes1.size();isp+=2) {
		std::pair<int,float> ind_pes1=pes1[isp];
		std::pair<int,float> ind_pes2=pes1[isp+1];
		float tot = ind_pes1.second+ind_pes2.second;
		if (tot>hitpe1)	  hitpe1=tot;
	      }
	      pes1CRT[j]=hitpe1;
	      hitpe1 = 0; 
	      for (uint isp=0;isp<pes2.size();isp+=2) {
		std::pair<int,float> ind_pes1=pes2[isp];
		std::pair<int,float> ind_pes2=pes2[isp+1];
		float tot = ind_pes1.second+ind_pes2.second;
		if (tot>hitpe1)	  hitpe1=tot;
	      }
	      pes2CRT[j]=hitpe1;	      
	    }

	    // Correct start and end points for space charge

	    geo::Point_t newStartP = startP; geo::Point_t newEndP = endP;
	    double xshift = vdrift*tzeroCRT[j];	  
	    if(sce->EnableCalSpatialSCE()) {

	      geo::Point_t fTrackPos = startP;  fTrackPos.SetX(startP.X()-xshift);
	      geo::Vector_t fPosOffsets = sce->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()});
	      newStartP = geo::Point_t{fTrackPos.X() - fPosOffsets.X(), fTrackPos.Y() + fPosOffsets.Y(), 
				       fTrackPos.Z() + fPosOffsets.Z()};
	      // std::cout << fPosOffsets.X() << " " <<   fPosOffsets.Y() << " " <<  fPosOffsets.Z() << std::endl;
	      
	      fTrackPos = endP;  fTrackPos.SetX(endP.X()-xshift);
	      fPosOffsets = sce->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()});
	      newEndP = geo::Point_t{fTrackPos.X() - fPosOffsets.X(), fTrackPos.Y() + fPosOffsets.Y(), 
				     fTrackPos.Z() + fPosOffsets.Z()};
	      
	      trkstartx[j]=newStartP.X();
	      trkstarty[j]=newStartP.Y();
	      trkstartz[j]=newStartP.Z();
	      trkendx[j]=newEndP.X();
	      trkendy[j]=newEndP.Y();
	      trkendz[j]=newEndP.Z();
	      // std::cout << fPosOffsets.X() << " " <<   fPosOffsets.Y() << " " <<  fPosOffsets.Z() << std::endl;
	    }
	    else {
	      newStartP.SetX(startP.X()-xshift); newEndP.SetX(endP.X()-xshift);
	    }
      
	    TVector3 trackstart(newStartP.X(),newStartP.Y(),newStartP.Z());
	    TVector3 trackend(newEndP.X(),newEndP.Y(),newEndP.Z());
	    
	    // angle between track extrap and CRT plane
	    TVector3 trackdir = trackstart-trackend;
	    TVector3 proj=trackdir;  
	    if (planeCRT[j]==0 || planeCRT[j]==3) proj.SetY(0.0);
	    else proj.SetX(0.0);
	    trackdir.SetMag(1.0); proj.SetMag(1.0);
	    double cosang=trackdir.Dot(proj);
	    if (cosang<0) cosang*=-1.0;
	    if (cosang>1) { if (verbose_) std::cout << "bad value cosang " << cosang << std::endl; cosang=0.0;}
	    angCRT[j]= acos(cosang);
	    
	  }}}// if there is a crt hit associated to the track
    }  // loop over tracks
  }   //  if (saveTPCtrackinfo)

  //get Optical Flash
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlash;
  evt.getByLabel(data_label_flash_, rawHandle_OpFlash);
  
  std::vector<recob::OpFlash> const& OpFlashCollection(*rawHandle_OpFlash);
  
  if(verbose_==1){ 
    std::cout<<"  OpFlashCollection.size()  "<<OpFlashCollection.size()<<std::endl; 
  }  //get Optical Flash
  

  //  fill tree
  run=frunNum;
  event=fEvtNum;
  subrun=fsubRunNum;
  

  // flashinfo
  nPMTflash = OpFlashCollection.size();
  if (nPMTflash>kMaxPMTflash) nPMTflash=kMaxPMTflash;
  for(int j = 0; j < nPMTflash; j++) {
    recob::OpFlash my_flash = OpFlashCollection[j];
    fl_time[j]=my_flash.Time();
    fl_pe[j]=my_flash.TotalPE();
    fl_y[j]=my_flash.YCenter();
    fl_z[j]=my_flash.ZCenter();
  }    
  

  //get CRTHits
  art::Handle< std::vector<crt::CRTHit> > rawHandle_hit;
  evt.getByLabel(data_labelhit_, rawHandle_hit); //
  
  //check to make sure the data we asked for is valid
  if(!rawHandle_hit.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTHits " << " in module " << data_labelhit_ << std::endl;
    std::cout << std::endl;
    return;
  }
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle_hit);
  if(verbose_==1){ 
    std::cout<<"  CRTHitCollection.size()  "<<CRTHitCollection.size()<<" label "<<data_labelhit_<<std::endl; 
    //  getchar();   
  }    //end get CRTHits


  nCRThits = CRTHitCollection.size();
  if (nCRThits>kMaxCRThits) nCRThits=kMaxCRThits;
  for(int j = 0; j < nCRThits; j++) {
    //std::cout<<"in CRT hits "<<j<<"\n";
    //fill tree
    crt::CRTHit my_CRTHit = CRTHitCollection[j];
    hit_time_s[j]=(double)my_CRTHit.ts0_s;
    hit_time0[j]=(double)my_CRTHit.ts0_ns - (double)evt_timeGPS_nsec + (double)fTimeZeroOffset;
    hit_time1[j]=(double)my_CRTHit.ts1_ns + (double)fHardDelay_; 
    hit_charge[j]=my_CRTHit.peshit;
    hit_plane[j]=(my_CRTHit.plane)%10;
    hit_posx[j]=my_CRTHit.x_pos;
    hit_posy[j]=my_CRTHit.y_pos;
    hit_posz[j]=my_CRTHit.z_pos;
    hit_feb1[j]=int(my_CRTHit.feb_id[0]);
    hit_feb2[j]=int(my_CRTHit.feb_id[1]);
    
    std::vector<std::pair<int,float>> pes = my_CRTHit.pesmap.find(hit_feb1[j])->second; 
    if (pes.size()==2) { // is simulated hit
      //std::cout<<"Sim CRT hits "<<j<<"\n";
      hit_mcflag[j]=1;
      //    std::map<uint8_t, std::vector<std::pair<int,float>>> pesmap;	      
      std::pair<int,float> ind_pes1 = pes[0];
      std::pair<int,float> ind_pes2 = pes[1];
      
      // std::cout << "first strip: feb  " << hit_feb1[j] << " map index feb " <<
      //   int(my_CRTHit.pesmap.find(hit_feb1[j])->first) << " sipm no " <<
      //   ind_pes1.first << " pes " << ind_pes1.second << " sipm no " <<
      //   ind_pes2.first << " pes " << ind_pes2.second << " strip no " <<
      //   int(0.5*(ind_pes1.first)) << std::endl;
      
      hit_strip1[j]=int(0.5*(ind_pes1.first));
      hit_pe1[j]=ind_pes1.second+ind_pes2.second;
      hit_sipm1a[j]=ind_pes1.second;    hit_sipm1b[j]=ind_pes2.second;
      pes = my_CRTHit.pesmap.find(hit_feb2[j])->second;
      ind_pes1 = pes[0]; ind_pes2 = pes[1];
      hit_strip2[j]=int(0.5*ind_pes1.first);
      hit_pe2[j]=ind_pes1.second+ind_pes2.second;
      hit_sipm2a[j]=ind_pes1.second;    hit_sipm2b[j]=ind_pes2.second;
      
      // if (hit_charge[j]>500.) hit_charge[j]=500.;
      // if (hit_pe1[j]>300.) hit_pe1[j]=300.;
      // if (hit_pe2[j]>300.) hit_pe2[j]=300.;
      // std::cout << "second strip: feb " << hit_feb2[j] << " map index feb " <<
      //   int(my_CRTHit.pesmap.find(hit_feb2[j])->first) << " sipm no " <<
      //   ind_pes1.first << " pes " << ind_pes1.second << " simp no " <<
      //   ind_pes2.first << " pes " << ind_pes2.second << " strip no " <<
      //   int(0.5*ind_pes1.first) << std::endl;
    }
    else {
      //  Data has info from all 32 sipms with no indication
      //   which ones were used to make this particular hit
      //
      // use the largest pe deposit on the feb, not necessarily correct
      //std::cout<<"Data CRT hits "<<j<<"\n";
      //    std::map<uint8_t, std::vector<std::pair<int,float>>> pesmap;	       //vector of pairs  
      hit_mcflag[j]=0;
      float hitpe1 = 0; int hitsipm1=-1;
      float spe1=0; float spe2=0;
      for (uint isp=0;isp<pes.size();isp+=2) {
	std::pair<int,float> ind_pes1=pes[isp];
	std::pair<int,float> ind_pes2=pes[isp+1];
	// std::cout << isp << " " << ind_pes1.first << " " << ind_pes1.second << std::endl;
	// std::cout << isp+1 << " " << ind_pes2.first << " " << ind_pes2.second << std::endl;
	float tot = ind_pes1.second+ind_pes2.second;
	if (tot>hitpe1) {
	  hitpe1=tot;
	  spe1=ind_pes1.second; spe2=ind_pes2.second;
	  hitsipm1=0.5*ind_pes1.first;
	}
      }
      hit_strip1[j]=-1;hit_pe1[j]=-1;
      if (hitsipm1>=0) {hit_strip1[j]=hitsipm1;      
	hit_sipm1a[j]=spe1;	hit_sipm1b[j]=spe2;
	hit_pe1[j]=hitpe1;}
      //
      std::vector<std::pair<int,float>> pes2 = my_CRTHit.pesmap.find(hit_feb2[j])->second; 
      float hitpe2 = 0; int hitsipm2=-1;
      spe1=0; spe2=0;
      for (uint isp=0;isp<pes2.size();isp+=2) {
	std::pair<int,float> ind_pes1=pes2[isp];
	std::pair<int,float> ind_pes2=pes2[isp+1];
	float tot = ind_pes1.second+ind_pes2.second;
	if (tot>hitpe2) {
	  hitpe2=tot;
	  spe1=ind_pes1.second; spe2=ind_pes2.second;
	  hitsipm2=0.5*ind_pes1.first;
	}
      }
      hit_strip2[j]=-1; hit_pe2[j]=-1;
      if (hitsipm2>=0) {hit_strip2[j]=hitsipm2;
	hit_sipm2a[j]=spe1;	hit_sipm2b[j]=spe2; 
	hit_pe2[j]=hitpe2;
      }
     //
	// std::cout << " feb1 strip1 pe1 " << hit_feb1[j] << " " << hit_strip1[j] << " " << 
	//   hit_pe1[j] << std::endl;
	// std::cout << " feb2 strip2 pe2 " << hit_feb2[j] << " " << hit_strip2[j] << " " << 
	//   hit_pe2[j] << std::endl;
      //
    } // if data
    

    //fillhistograms
    int thisplane = (my_CRTHit.plane%10);
    if (thisplane==0) HitDistBot->Fill(my_CRTHit.z_pos,my_CRTHit.x_pos);
    else if (thisplane==1) HitDistFT->Fill(my_CRTHit.z_pos,my_CRTHit.y_pos);
    else if (thisplane==2) HitDistPipe->Fill(my_CRTHit.z_pos,my_CRTHit.y_pos);
    else if (thisplane==3) HitDistTop->Fill(my_CRTHit.z_pos,my_CRTHit.x_pos);
  }// loop over hits


  //get CRT tzeros
  art::Handle< std::vector<crt::CRTTzero> > rawHandle_tz;
  evt.getByLabel(data_labeltzero_, rawHandle_tz); 
  
  //check to make sure the data we asked for is valid
  if(!rawHandle_tz.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTtzeros " << " in module " << data_labeltzero_ << std::endl;
    std::cout << std::endl;
    return;
  }
  std::vector<crt::CRTTzero> const& CRTtzCollection(*rawHandle_tz);
  if(verbose_==1){ 
    std::cout<<"  CRTtzCollection.size()  "<<CRTtzCollection.size()<<std::endl; 
    //  getchar();   
  }    //end get CRTTzeros


  nCRTtzeros = CRTtzCollection.size();
  if (nCRTtzeros>kMaxCRTtzs) nCRTtzeros=kMaxCRTtzs;
  for(int j = 0; j < nCRTtzeros; j++) {
    crt::CRTTzero my_CRTtz = CRTtzCollection[j];
    tz_time_s[j]=(double)my_CRTtz.ts0_s;
    tz_time0[j]=(double)my_CRTtz.ts0_ns - (double)evt_timeGPS_nsec + (double)fTimeZeroOffset;
    tz_time1[j]=(double)my_CRTtz.ts1_ns + (double)fHardDelay_; 
    tz_hits0[j]=my_CRTtz.nhits[0];
    tz_hits1[j]=my_CRTtz.nhits[1];
    tz_hits2[j]=my_CRTtz.nhits[2];
    tz_hits3[j]=my_CRTtz.nhits[3];
    tz_pes0[j]=my_CRTtz.pes[0];
    tz_pes1[j]=my_CRTtz.pes[1];
    tz_pes2[j]=my_CRTtz.pes[2];
    tz_pes3[j]=my_CRTtz.pes[3];
  }


  //get CRTTracks
  art::Handle< std::vector<crt::CRTTrack> > rawHandle_track;
  evt.getByLabel(data_labeltrack_, rawHandle_track); 
  
  //check to make sure the data we asked for is valid
  if(!rawHandle_track.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTTracks " << " in module " << data_labeltrack_ << std::endl;
    std::cout << std::endl;
    return;
  }

  std::vector<crt::CRTTrack> const& CRTTrackCollection(*rawHandle_track);
  if(verbose_==1){ 
    std::cout<<"  CRTTrackCollection.size()  "<<CRTTrackCollection.size()<<std::endl; 
    getchar();   
  }


  nCRTtracks = CRTTrackCollection.size();
  if (nCRTtracks>kMaxCRTtracks) nCRTtracks=kMaxCRTtracks;

  for(int j = 0; j <nCRTtracks; j++) {
    crt::CRTTrack my_CRTTrack = CRTTrackCollection[j];
    double temp = (my_CRTTrack.y1_pos-my_CRTTrack.y2_pos)*(my_CRTTrack.y1_pos-my_CRTTrack.y2_pos)
      +(my_CRTTrack.x1_pos-my_CRTTrack.x2_pos)*(my_CRTTrack.x1_pos-my_CRTTrack.x2_pos);
    double thetatemp =  atan2(sqrt(temp),my_CRTTrack.z1_pos-my_CRTTrack.z2_pos);
    double phitemp = atan2(my_CRTTrack.y1_pos-my_CRTTrack.y2_pos,my_CRTTrack.x1_pos-my_CRTTrack.x2_pos);
    // flip track to point downwards if needed
    if (phitemp<0) { ct_phi[j]=phitemp+3.14159; ct_theta[j]=3.14159-thetatemp;}
    else { ct_theta[j]=thetatemp;     ct_phi[j]=phitemp;}
    ct_length[j]=my_CRTTrack.length;
    ct_time_sec[j]=(double)my_CRTTrack.ts0_s;
    ct_time0[j]=(double)my_CRTTrack.ts0_ns;
    ct_time1[j]=(double)my_CRTTrack.ts1_ns;
    // ct_time0[j]=(double)my_CRTTrack.ts0_ns - (double)evt_timeGPS_nsec;
    // ct_time1[j]=(double)my_CRTTrack.ts1_ns + (double)fHardDelay_;    + 40000 for hardware offset;;
    // std::cout << my_CRTTrack.ts0_ns << " " << ct_time0[j] << std::endl;
    // std::cout << my_CRTTrack.ts1_ns << " " << ct_time1[j] << std::endl;
    ct_x1[j]=my_CRTTrack.x1_pos;
    ct_y1[j]=my_CRTTrack.y1_pos;
    ct_z1[j]=my_CRTTrack.z1_pos;
    ct_x2[j]=my_CRTTrack.x2_pos;
    ct_y2[j]=my_CRTTrack.y2_pos;
    ct_z2[j]=my_CRTTrack.z2_pos;

    //fill histograms
    hplavspla->Fill((my_CRTTrack.plane1)%10,(my_CRTTrack.plane2)%10);
    hTlength->Fill(my_CRTTrack.length);
    double time_diff = my_CRTTrack.ts0_ns_h1-my_CRTTrack.ts0_ns_h2;
    double time_diffABS = fabs(time_diff);
    hTtime->Fill(time_diffABS);
    hTlengthvsTimeAbs->Fill(my_CRTTrack.length,time_diffABS);
    hTlengthvsTimeAbs_prof->Fill(my_CRTTrack.length,time_diffABS);
    hTlengthvsTime->Fill(my_CRTTrack.length,time_diff);
    htheta->Fill(57.30*my_CRTTrack.thetaxy);
    if (my_CRTTrack.phizy>3.14159) 
      hphi->Fill(57.30*(my_CRTTrack.phizy-3.14159));
    else    hphi->Fill(57.30*my_CRTTrack.phizy);
    hts0_ns->Fill(my_CRTTrack.ts0_ns);


  }


  fTree->Fill();

  
}

void TrackDump::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("trackdump","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("nCRThits",&nCRThits,"nCRThits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nCRThits]/I");
  fTree->Branch("hit_time_s",hit_time_s,"hit_time_s[nCRThits]/D");
  fTree->Branch("hit_time0",hit_time0,"hit_time0[nCRThits]/D");
  fTree->Branch("hit_time1",hit_time1,"hit_time1[nCRThits]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nCRThits]/D");
  fTree->Branch("hit_posx",hit_posx,"hit_posx[nCRThits]/D");
  fTree->Branch("hit_posy",hit_posy,"hit_posy[nCRThits]/D");
  fTree->Branch("hit_posz",hit_posz,"hit_posz[nCRThits]/D");
  fTree->Branch("hit_feb1",hit_feb1,"hit_feb1[nCRThits]/I");
  fTree->Branch("hit_feb2",hit_feb2,"hit_feb2[nCRThits]/I");
  fTree->Branch("hit_strip1",hit_strip1,"hit_strip1[nCRThits]/I");
  fTree->Branch("hit_strip2",hit_strip2,"hit_strip2[nCRThits]/I");
  fTree->Branch("hit_pe1",hit_pe1,"hit_pe1[nCRThits]/D");
  fTree->Branch("hit_pe2",hit_pe2,"hit_pe2[nCRThits]/D");
  fTree->Branch("hit_sipm1a",hit_sipm1a,"hit_sipm1a[nCRThits]/D");
  fTree->Branch("hit_sipm1b",hit_sipm1b,"hit_sipm1b[nCRThits]/D");
  fTree->Branch("hit_sipm2a",hit_sipm2a,"hit_sipm2a[nCRThits]/D");
  fTree->Branch("hit_sipm2b",hit_sipm2b,"hit_sipm2b[nCRThits]/D");
  fTree->Branch("hit_mcflag",hit_mcflag,"hit_mcflag[nCRThits]/I");
  // CRT tzeros
  fTree->Branch("nCRTtzeros",&nCRTtzeros,"nCRTtzeros/I");
  fTree->Branch("tz_time_s",tz_time_s,"tz_time_s[nCRTtzeros]/D");
  fTree->Branch("tz_time0",tz_time0,"tz_time0[nCRTtzeros]/D");
  fTree->Branch("tz_time1",tz_time1,"tz_time1[nCRTtzeros]/D");
  fTree->Branch("tz_hits0",tz_hits0,"tz_hits0[nCRTtzeros]/I");
  fTree->Branch("tz_hits1",tz_hits1,"tz_hits1[nCRTtzeros]/I");
  fTree->Branch("tz_hits2",tz_hits2,"tz_hits2[nCRTtzeros]/I");
  fTree->Branch("tz_hits3",tz_hits3,"tz_hits3[nCRTtzeros]/I");
  fTree->Branch("tz_pes0",tz_pes0,"tz_pes0[nCRTtzeros]/D");
  fTree->Branch("tz_pes1",tz_pes1,"tz_pes1[nCRTtzeros]/D");
  fTree->Branch("tz_pes2",tz_pes2,"tz_pes2[nCRTtzeros]/D");
  fTree->Branch("tz_pes3",tz_pes3,"tz_pes3[nCRTtzeros]/D");
  // CRT tracks
  fTree->Branch("nCRTtracks",&nCRTtracks,"nCRTtracks/I");
  fTree->Branch("ct_theta",ct_theta,"ct_theta[nCRTtracks]/D");
  fTree->Branch("ct_phi",ct_phi,"ct_phi[nCRTtracks]/D");
  fTree->Branch("ct_length",ct_length,"ct_length[nCRTtracks]/D");
  fTree->Branch("ct_time_sec",ct_time_sec,"ct_time_sec[nCRTtracks]/D");
  fTree->Branch("ct_time0",ct_time0,"ct_time0[nCRTtracks]/D");
  fTree->Branch("ct_time1",ct_time1,"ct_time1[nCRTtracks]/D");
  fTree->Branch("ct_x1",ct_x1,"ct_x1[nCRTtracks]/D");
  fTree->Branch("ct_y1",ct_y1,"ct_y1[nCRTtracks]/D");
  fTree->Branch("ct_z1",ct_z1,"ct_z1[nCRTtracks]/D");
  fTree->Branch("ct_x2",ct_x2,"ct_x2[nCRTtracks]/D");
  fTree->Branch("ct_y2",ct_y2,"ct_y2[nCRTtracks]/D");
  fTree->Branch("ct_z2",ct_z2,"ct_z2[nCRTtracks]/D");
  //TPC tracks
  if (fSaveTPCTrackInfo) {
  fTree->Branch("nTPCtracks",&nTPCtracks,"nTPCtracks/I");
  fTree->Branch("trkstartx",trkstartx,"trkstartx[nTPCtracks]/D");
  fTree->Branch("trkstarty",trkstarty,"trkstarty[nTPCtracks]/D");
  fTree->Branch("trkstartz",trkstartz,"trkstartz[nTPCtracks]/D");
  fTree->Branch("trkendx",trkendx,"trkendx[nTPCtracks]/D");
  fTree->Branch("trkendy",trkendy,"trkendy[nTPCtracks]/D");
  fTree->Branch("trkendz",trkendz,"trkendz[nTPCtracks]/D");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[nTPCtracks]/D");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[nTPCtracks]/D");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[nTPCtracks]/D");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[nTPCtracks]/D");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[nTPCtracks]/D");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[nTPCtracks]/D");
  fTree->Branch("trktheta",trktheta,"trktheta[nTPCtracks]/D");
  fTree->Branch("trkphi",trkphi,"trkphi[nTPCtracks]/D");
  fTree->Branch("trklen",trklen,"trklen[nTPCtracks]/D");
  fTree->Branch("tzeroACPT",tzeroACPT,"tzeroACPT[nTPCtracks]/D");
  fTree->Branch("tzeroCRT",tzeroCRT,"tzeroCRT[nTPCtracks]/D");
  fTree->Branch("planeCRT",planeCRT,"planeCRT[nTPCtracks]/I");
  fTree->Branch("feb1CRT",feb1CRT,"feb1CRT[nTPCtracks]/I");
  fTree->Branch("feb2CRT",feb2CRT,"feb2CRT[nTPCtracks]/I");
  fTree->Branch("pes1CRT",pes1CRT,"pes1CRT[nTPCtracks]/D");
  fTree->Branch("pes2CRT",pes2CRT,"pes2CRT[nTPCtracks]/D");
  fTree->Branch("xCRT",xCRT,"xCRT[nTPCtracks]/D");
  fTree->Branch("yCRT",yCRT,"yCRT[nTPCtracks]/D");
  fTree->Branch("zCRT",zCRT,"zCRT[nTPCtracks]/D");
  fTree->Branch("angCRT",angCRT,"angCRT[nTPCtracks]/D");
  }
  if (fSavePMTFlashInfo) {
    fTree->Branch("nPMTflash",&nPMTflash,"nPMTflash/I");
    fTree->Branch("fl_time",fl_time,"fl_time[nPMTflash]/D");
    fTree->Branch("fl_pe",fl_pe,"fl_pe[nPMTflash]/D");
    fTree->Branch("fl_y",fl_y,"fl_y[nPMTflash]/D");
    fTree->Branch("fl_z",fl_z,"fl_z[nPMTflash]/D");
  }

  hplavspla = tfs->make<TH2F>("hplavspla","PlanevsPlane",4,0,4,4,0,4);
  hplavspla->GetXaxis()->SetTitle("Plane (0=Bottom, 1=FT, 2=Pipe, 3=Top)");
  hplavspla->GetYaxis()->SetTitle("Plane (0=Bottom, 1=FT, 2=Pipe, 3=Top)");
  hplavspla->GetZaxis()->SetTitle("Entries/bin");
  hplavspla->SetOption("COLZ");

  hTvsH = tfs->make<TH2F>("hTvsH","Track_vs_Hits",500,0,500,500,0,500);
  hTvsH->GetXaxis()->SetTitle("Number of CRTHits per event");
  hTvsH->GetYaxis()->SetTitle("Number of CRTTracks per event");
  hTvsH->GetZaxis()->SetTitle("Entries/bin");
  hTvsH->SetOption("COLZ");

  hTlength = tfs->make<TH1F>("hTlength","Track_Length",1500,0,1500);
  hTlength->GetXaxis()->SetTitle("Track_Length (cm)");
  hTlength->GetYaxis()->SetTitle("Entries/bin");

  hTtime = tfs->make<TH1F>("hTtime","Track_time",120,-10,110);
  hTtime->GetXaxis()->SetTitle("Track_time (ns)");
  hTtime->GetYaxis()->SetTitle("Entries/bin");

  hTlengthvsTime = tfs->make<TH2F>("hTlengthvsTime","Track_LengthvsTime",1500,0,1500,200,-100,100);
  hTlengthvsTime->GetXaxis()->SetTitle("Track_Length (cm)");
  hTlengthvsTime->GetYaxis()->SetTitle("Track_time (ns)");
  hTlengthvsTime->GetZaxis()->SetTitle("Entries/bin");
  hTlengthvsTime->SetOption("COLZ");

  hTlengthvsTimeAbs = tfs->make<TH2F>("hTlengthvsTimeAbs","Track_LengthvsTimeAbs",1500,0,1500,110,-10,100);
  hTlengthvsTimeAbs->GetXaxis()->SetTitle("Track_Length (cm)");
  hTlengthvsTimeAbs->GetYaxis()->SetTitle("Track_time (ns)");
  hTlengthvsTimeAbs->GetZaxis()->SetTitle("Entries/bin");
  hTlengthvsTimeAbs->SetOption("COLZ");

  hTlengthvsTimeAbs_prof = tfs->make<TProfile>("hTlengthvsTimeAbs_prof","Track_LengthvsTimeAbs_prof",1500,0,1500,"s");
  hTlengthvsTimeAbs_prof->GetXaxis()->SetTitle("Track_Length (cm)");
  hTlengthvsTimeAbs_prof->GetYaxis()->SetTitle("Track_time (ns)");

  htheta = tfs->make<TH1F>("htheta","Track_theta",900,0,180);
  htheta->GetXaxis()->SetTitle("Theta_xy (º)");
  htheta->GetYaxis()->SetTitle("Entries/bin");
 
  hphi = tfs->make<TH1F>("hphi","Track_phi",900,0,180);
  hphi->GetXaxis()->SetTitle("Phi_zy (º)");
  hphi->GetYaxis()->SetTitle("Entries/bin");

  hts0_ns = tfs->make<TH1F>("hts0_ns","Track_time_ns",100000,0,1e9);
  hts0_ns->GetXaxis()->SetTitle("Track time (ns)");
  hts0_ns->GetYaxis()->SetTitle("Entries/bin");


  double inch =2.54; //inch in cm                                                                                                                          
  HitDistBot = tfs->make<TH2F>("hBottom","Bottom",125,-700+205*inch,-700+205*inch+125*10.89,60,-300+50.4*inch,-300+50.4*inch+60*10.89);
  HitDistBot->GetXaxis()->SetTitle("Length long the beam (cm)");
  HitDistBot->GetYaxis()->SetTitle("Length along the drift (cm)");
  HitDistBot->GetZaxis()->SetTitle("Entries/bin");
  HitDistBot->SetOption("COLZ");

  HitDistFT = tfs->make<TH2F>("hFeedthroughSide","Feedthrough Side",125,-704+205*inch,-704+205*inch+125*10.89,60,-308-19.1*inch,-308-19.1*inch+60*10.89);
  HitDistFT->GetXaxis()->SetTitle("Length along the beam (cm)");
  HitDistFT->GetYaxis()->SetTitle("Height (cm)");
  HitDistFT->GetZaxis()->SetTitle("Entries/bin");
  HitDistFT->SetOption("COLZ");

  HitDistPipe = tfs->make<TH2F>("hPipeSide","Pipe Side",125,-704+205*inch,-704+205*inch+125*10.89,60,-294-19.1*inch,-294-19.1*inch+60*10.89);
  HitDistPipe->GetXaxis()->SetTitle("Length along the beam (cm)");
  HitDistPipe->GetYaxis()->SetTitle("Height (cm)");
  HitDistPipe->GetZaxis()->SetTitle("Entries/bin");
  HitDistPipe->SetOption("COLZ");

  HitDistTop = tfs->make<TH2F>("hTop","Top",125,-701+205*inch,-701+205*inch+125*11.38,80,2-170-300+50.4*inch,2-170-300+50.4*inch+80*11.38);
  HitDistTop->GetXaxis()->SetTitle("Length along the beam (cm)");
  HitDistTop->GetYaxis()->SetTitle("Length along the drift (cm)");
  HitDistTop->GetZaxis()->SetTitle("Entries/bin");
  HitDistTop->SetOption("COLZ");



}

void TrackDump::endJob()
{
  // Implementation of optional member function here.
  
  
  //fTree->Write();
  
}


void TrackDump::ResetVars()
{
  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  nCRThits = 0;
  for (int i = 0; i<kMaxCRThits; ++i){
    hit_plane[i] = -999;
    hit_time_s[i] = -99999.;
    hit_time0[i] = -99999.;
    hit_time1[i] = -99999.;
    hit_charge[i] = -99999.;
    hit_posx[i] = -99999.;
    hit_posy[i] = -99999.;
    hit_posz[i] = -99999.;
    hit_feb1[i] = -999;
    hit_feb2[i] = -999;
    hit_strip1[i] = -1;
    hit_strip2[i] = -1;
    hit_pe1[i] = -999;
    hit_pe2[i] = -999;
    hit_mcflag[i] = -1;
  }


  nCRTtracks=0;
  for (int j = 0; j<kMaxCRTtracks; ++j){
    ct_theta[j]=-99999.;
    ct_phi[j]=-99999.;
    ct_length[j]=-99999.;
    ct_time_sec[j]=-99999.;
    ct_time0[j]=-99999.;
    ct_time1[j]=-99999.;
    ct_x1[j]=-99999.;
    ct_y1[j]=-99999.;
    ct_z1[j]=-99999.;
    ct_x2[j]=-99999.;
    ct_y2[j]=-99999.;
    ct_z2[j]=-99999.;
  }

  if (fSaveTPCTrackInfo) {
  nTPCtracks=0;
  for (int i = 0; i<kMaxTPCtracks; ++i){
    trkstartx[i]=-9999.;
    trkstarty[i]=-9999.;
    trkstartz[i]=-9999.;
    trkendx[i]=-9999.;
    trkendy[i]=-9999.;
    trkendz[i]=-9999.;
    trkstartdcosx[i]=-9999.;
    trkstartdcosy[i]=-9999.;
    trkstartdcosz[i]=-9999.;
    trkenddcosx[i]=-9999.;
    trkenddcosy[i]=-9999.;
    trkenddcosz[i]=-9999.;
    trktheta[i]=-9999.;
    trkphi[i]=-9999.;
    trklen[i]=-9999.;
    tzeroACPT[i]=-9999.;
    tzeroCRT[i]=-9999.;
  //  double dcaCRT[kMaxTPCtracks];
    planeCRT[i]=-1;
    feb1CRT[i]=-1;
    feb2CRT[i]=-1;
    pes1CRT[i]=-9999.;
    pes2CRT[i]=-9999.;
    xCRT[i]=-9999.;
    yCRT[i]=-9999.;
    zCRT[i]=-9999.;
    angCRT[i]=-9999.;
  }
}

}
DEFINE_ART_MODULE(TrackDump)


