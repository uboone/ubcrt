////////////////////////////////////////////////////////////////////////
// Class:       CRTHitPositionPerModule
// Module Type: analyzer
// File:        CRTHitPositionPerModule_module.cc
//
// Generated at Mon Sept  10 13:49:00 2018 by Elena Gramellini
// Scope of this analyszer is dumping out all the info needed to 
// Compute the rate of cosmic rays in the CRT top module
// 
// [   ] Read out CRT hits
// [   ] Identify Corresponding Module
// [   ] Is this data or MC?
// [   ] Delta T consecutive hits in same module
// [   ] Hit X Position 
// [   ] Hit Y Position 
// [   ] Hit Z Position 
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

#include "art_root_io/TFileService.h"

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTTzero.hh"
#include "ubobj/CRT/CRTTrack.hh"
#include "ubcrt/CRT/CRTAuxFunctions.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

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
//const int kMaxCRTtzeros = 1000;
//const int kMaxCRTtracks = 1000;
//const int kMaxTPCtracks = 100;
//const int kMaxPMTflashes = 100;


 // namespace crt {
 //   class CRTHitPositionPerModule;
 // }

class CRTHitPositionPerModule : public art::EDAnalyzer {
public:
  explicit CRTHitPositionPerModule(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTHitPositionPerModule(CRTHitPositionPerModule const &) = delete;
  CRTHitPositionPerModule(CRTHitPositionPerModule &&) = delete;
  CRTHitPositionPerModule & operator = (CRTHitPositionPerModule const &) = delete;
  CRTHitPositionPerModule & operator = (CRTHitPositionPerModule &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  art::ServiceHandle<art::TFileService> tfs;

  // Declare member data here.
  double nEvt = 0 ;
  //  uint32_t fEvtNum; //Number of current event                       
  //uint32_t frunNum;                //Run Number taken from event  
  //uint32_t fsubRunNum;             //Subrun Number taken from event         
  std::string  fTrackModuleLabel;
  bool fSaveTPCTrackInfo;
  bool fSavePMTFlashInfo;
  std::string data_labeltrack_;
  std::string data_labeltzero_;
  std::string data_labelhit_;
  std::string data_label_flash_;
  std::string data_label_DAQHeader_;
  int fHardDelay_;
  int verbose_;
  bool isData;


  double evttime_sec;
  double evttime_nsec;
  double evttime_GPS_sec;
  double evttime_GPS_nsec;
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
  // CRT tzeros
  //  int nCRTtzeros;
  //double tz_time_s[kMaxCRTtzeros];
  //double tz_time0[kMaxCRTtzeros];
  //double tz_time1[kMaxCRTtzeros];

  
  TH2F* HitDistBot;
  TH2F* HitDistFT;
  TH2F* HitDistPipe;
  TH2F* HitDistTop;
  
  TH1F* hHitX_Module[73];
  TH1F* hHitY_Module[73];
  TH1F* hHitZ_Module[73];



};


CRTHitPositionPerModule::CRTHitPositionPerModule(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    fTrackModuleLabel(p.get<std::string>("TrackModuleLabel")),
    fSaveTPCTrackInfo(p.get< bool >("SaveTPCTrackInfo", false)), 
    fSavePMTFlashInfo(p.get< bool >("SavePMTFlashInfo", false)), 
    data_labeltrack_(p.get<std::string>("data_labeltrack")),
    data_labeltzero_(p.get<std::string>("data_labeltzero")),
    data_labelhit_(p.get<std::string>("data_labelhit")),
    data_label_flash_(p.get<std::string>("data_label_flash_")),
    data_label_DAQHeader_(p.get<std::string>("data_label_DAQHeader_")),
    fHardDelay_(p.get<int>("fHardDelay",40000)), // 40 us
    verbose_(p.get<int>("verbose")),
    isData(p.get<bool>("isData",false))
    // More initializers here.    
{
}

void CRTHitPositionPerModule::analyze(art::Event const & evt)
{
  nEvt++;
  
  art::Timestamp        evtTime  = evt.time();
  auto             evt_time_sec  = evtTime.timeHigh();
  auto             evt_time_nsec = evtTime.timeLow();

  double         evt_timeGPS_sec   = 0;
  double         evt_timeGPS_nsec  = 0;
  
  if (isData)
    {
      //get DAQ Header                                                                  
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
      art::Timestamp evtTimeGPS        = my_DAQHeader.gps_time();          // This is the GPS time as taken from the DAQ  
      evt_timeGPS_sec   = evtTimeGPS.timeHigh(); 
      evt_timeGPS_nsec  = (double)evtTimeGPS.timeLow();
      //      art::Timestamp evtTimeNTP        = my_DAQHeader.ntp_time();          // This is the NTP time as taken from the DAQ  
    }//isData?

  
  //fill tree variables
  evttime_sec      = evt_time_sec;
  evttime_nsec     = evt_time_nsec;
  evttime_GPS_sec  = evt_timeGPS_sec;
  evttime_GPS_nsec = evt_timeGPS_nsec;

  
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

  nCRThits = CRTHitCollection.size();
  if (nCRThits>kMaxCRThits) nCRThits=kMaxCRThits;
  for(int j = 0; j < nCRThits; j++) {
    //fill tree
    crt::CRTHit my_CRTHit = CRTHitCollection[j];
    hit_time_s[j]         = (double)my_CRTHit.ts0_s;
    hit_time0[j]          = (double)my_CRTHit.ts0_ns - (double)evt_timeGPS_nsec; // Time in the TPC reference frame (or PMT) 
    hit_time1[j]          = (double)my_CRTHit.ts1_ns + (double)fHardDelay_;      // Time in the TPC reference frame (or PMT)
    hit_charge[j]         = my_CRTHit.peshit;
    hit_plane[j]          = my_CRTHit.plane;
    hit_posx[j]           = my_CRTHit.x_pos;
    hit_posy[j]           = my_CRTHit.y_pos;
    hit_posz[j]           = my_CRTHit.z_pos;
    
    // Bottom   
    if (my_CRTHit.plane == 0) HitDistBot->Fill(my_CRTHit.z_pos,my_CRTHit.x_pos); 
    // FT Side
    if (my_CRTHit.plane == 1) HitDistFT->Fill(my_CRTHit.z_pos,my_CRTHit.y_pos);
    // Pipe Side
    if (my_CRTHit.plane == 2) HitDistPipe->Fill(my_CRTHit.z_pos,my_CRTHit.y_pos);
    // Top
    if (my_CRTHit.plane == 3) HitDistTop->Fill(my_CRTHit.z_pos,my_CRTHit.x_pos);
    
    // Every hit is composed by 2 FEBs: the horizontal and vertical one
    auto febID_v = my_CRTHit.feb_id;
    // Loop on the vector that composed the hit
    for (auto const febID : febID_v) 
      {  
	// Let's find the right histograms
	for (size_t iHist = 0; iHist < 73; iHist++)
	  {
	    // Check if title contains FEB STRING
	    // 1) convert uint_8 into string
	    std::ostringstream convert;
	    convert << (int) febID;
	    std::string key_string = convert.str();
	    // 2) get the histo title and chop it up 
	    std::string title = hHitX_Module[iHist]->GetTitle();
	    std::istringstream ss(title);
	    std::string token;
	    bool isThisTheHisto   = false;
	    while(std::getline(ss, token, '_'))
	      // 3) check if the feb name is contained and the number of characters is correct
	      if (token.find(key_string) !=std::string::npos  && key_string.size() == token.size())
		isThisTheHisto = true;
	    
	    if (isThisTheHisto)
	      {
		hHitX_Module[iHist] ->Fill(hit_posx[j]);
		hHitY_Module[iHist] ->Fill(hit_posy[j]);
		hHitZ_Module[iHist] ->Fill(hit_posz[j]);
	      }// Am I in the right histogram?
	    
	    isThisTheHisto = false;
	  }
      }// Loop on FEB vector
    
  }//Loop on CRT hits
 
}

void CRTHitPositionPerModule::beginJob()
{
  // Implementation of optional member function here.
  
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
  HitDistTop->GetXaxis()->SetTitle("Lenght along the beam (cm)");
  HitDistTop->GetYaxis()->SetTitle("Lenght along the drift (cm)");
  HitDistTop->GetZaxis()->SetTitle("Entries/bin");
  HitDistTop->SetOption("COLZ");

  
  
    
  std::string febNames[73] = {"11"  , "12" ,"14" ,"17"  ,"18"  ,"19"  ,  "22" ,"23"  ,  "24",
			      "105" ,"106" ,"107","108" ,"109" ,"111" , "112" ,"113" , "114", 
			      "115" ,"116" ,"117","118" ,"119" ,"120" , "121" ,"123" , "124", 
			      "125" ,"126" ,"127","128" ,"129" ,"195" ,  "26" , "27" , "28" , 
			      "29"  , "30" , "31", "52" , "56" , "57" , "58"  , "59" , "60" , 
			      "61"  , "15" , "16", "20" , "21" , "32" , "33"  , "34" , "35" , 
			      "36"  , "37" , "38", "39" , "40" , "41" , "42" ,  "43" , "44" , 
			      "45"  , "46" , "47", "48" , "49" , "50" , "51" ,  "53" , "54" , "55"};

  for (size_t i = 0; i < 73; i++)
    { 
      //Just naming this stuff... It's a nightmare
      // X Position
      std::string  hHitX_Module_str1  = "hHitX_Module_"  + febNames[i];
      std::string  hHitX_Module_str2  = "hHitX_Module_"  + febNames[i] + "; Hit X [cm]; Count";
      const char * hHitX_Module_Name1 = new char(hHitX_Module_str1.size()); 
      const char * hHitX_Module_Name2 = new char(hHitX_Module_str2.size()); 
      hHitX_Module_Name1 = hHitX_Module_str1.c_str();
      hHitX_Module_Name2 = hHitX_Module_str2.c_str();
      hHitX_Module[i] = tfs->make<TH1F>(hHitX_Module_Name1 ,hHitX_Module_Name2,  800,-300,500);

      // Y Position
      std::string  hHitY_Module_str1  = "hHitY_Module_"  + febNames[i];
      std::string  hHitY_Module_str2  = "hHitY_Module_"  + febNames[i] + "; Hit Y [cm]; Count";
      const char * hHitY_Module_Name1 = new char(hHitY_Module_str1.size()); 
      const char * hHitY_Module_Name2 = new char(hHitY_Module_str2.size()); 
      hHitY_Module_Name1 = hHitY_Module_str1.c_str();
      hHitY_Module_Name2 = hHitY_Module_str2.c_str();
      hHitY_Module[i] = tfs->make<TH1F>(hHitY_Module_Name1 ,hHitY_Module_Name2,  1000,-300,700);

      // Z Position
      std::string  hHitZ_Module_str1  = "hHitZ_Module_"  + febNames[i];
      std::string  hHitZ_Module_str2  = "hHitZ_Module_"  + febNames[i] + "; Hit Z [cm]; Count";
      const char * hHitZ_Module_Name1 = new char(hHitZ_Module_str1.size()); 
      const char * hHitZ_Module_Name2 = new char(hHitZ_Module_str2.size()); 
      hHitZ_Module_Name1 = hHitZ_Module_str1.c_str();
      hHitZ_Module_Name2 = hHitZ_Module_str2.c_str();
      hHitZ_Module[i] = tfs->make<TH1F>(hHitZ_Module_Name1 ,hHitZ_Module_Name2, 1500,-300,1200);
    }

 

}

void CRTHitPositionPerModule::endJob()
{  
  std::cout<<"N events "<<nEvt<<"\n";
}


DEFINE_ART_MODULE(CRTHitPositionPerModule)


