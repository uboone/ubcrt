////////////////////////////////////////////////////////////////////////
// Class:       CRTDataQuality
// Module Type: analyzer
// File:        CRTDataQuality_module.cc
//
// Generated at Thur March 28 2019 by Elena Gramellini
// Scope of this analyzer is a simple data quality monitor
// Compute the rate of cosmic rays in all CRT modules per date
// 
// [ x ] Read out CRT hits
// [ x ] Identify Corresponding Module
// [ x ] Count N Hit Per module in DeltaT
// [   ] Store the following in a ttree:
//            [ x ] date
//            [ x ] N hit per module
//            [ x ] FEBIndex
//            [ x ] DeltaT Readout
//            [ x ] AvgPe
//            [   ] AvgX, AvgY, AvgZ count
//            [   ] MaxX, MaxY, MaxZ count
//            [   ] MaxX, MaxY, MaxZ position
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
#include <iomanip>
#include <math.h>
#include <time.h>       /* time_t, struct tm, time, localtime, strftime */



//const int kMaxCRThits = 1000;
//const int kMaxCRTtzeros = 1000;
//const int kMaxCRTtracks = 1000;
//const int kMaxTPCtracks = 100;
//const int kMaxPMTflashes = 100;


 // namespace crt {
 //   class CRTDataQuality;
 // }

class CRTDataQuality : public art::EDAnalyzer {
public:
  explicit CRTDataQuality(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTDataQuality(CRTDataQuality const &) = delete;
  CRTDataQuality(CRTDataQuality &&) = delete;
  CRTDataQuality & operator = (CRTDataQuality const &) = delete;
  CRTDataQuality & operator = (CRTDataQuality &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void ResetVar();
  void calculateFlashMatch(std::vector<crt::CRTHit> const CRTHitCollection,    std::vector<recob::OpFlash> const OpFlashCollection, double event_timeGPS_ns, 
			   int &nFlashes, int &nMatchedFlashes, std::vector<double> &CRTFlashMatchTimeDiff);
  void calculateTrackMatch(   art::Handle< std::vector<recob::Track>  > trackListHandle,  art::FindMany<anab::T0> trk_t0C_assn_v, int Tcollsize, int &nTracks , int &nMatchedTracks );

private:

  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<art::TFileService> tfs1;

  // Declare member data here.
  std::string data_labelhit_;
  std::string data_label_DAQHeader_;
  std::string data_label_TPCtrack_;
  std::string data_label_flash_;
  std::string data_label_match_;
  double minT_ ;
  double maxT_ ;
  double flash_match_timecut;  // in us
  double min_length_TPCtrack;
  int    fTimeZeroOffset; 
  bool   applyTimeOffSet_;
  bool   verbose_;

  // TTree Variables
  TH1D*  hHitTime;
  TH1D*  hHitTimeAfterCut;
  TH1D*  hModule_X[73];
  TH1D*  hModule_Y[73];
  TH1D*  hModule_Z[73];
  TH1D* hFlashDiff;
  TH1D* hTrackDCA;
  TH1D* hTrackDCA3;

  TTree* fTree;
  int run;
  int subrun;
  int event;
  int date; // Time in seconds from linux start time
  int nFlashes;
  int nMatchedFlashes;
  int nTracks;
  int nMatchedTracks;
  std::vector<double> CRTFlashMatchTimeDiff;
  // CRT Modules
  int nCRThits[73];
  double AvgPe[73];


  double reaoutTime;
  int febIndex[73]  = {11  , 12 ,14 ,17  ,18  ,19  ,  22 ,23  ,  24,
		       105 ,106 ,107,108 ,109 ,111 , 112 ,113 , 114, 
		       115 ,116 ,117,118 ,119 ,120 , 121 ,123 , 124, 
		       125 ,126 ,127,128 ,129 ,195 ,  26 , 27 , 28 , 
		       29  , 30 , 31, 52 , 56 , 57 , 58  , 59 , 60 , 
		       61  , 15 , 16, 20 , 21 , 32 , 33  , 34 , 35 , 
		       36  , 37 , 38, 39 , 40 , 41 , 42 ,  43 , 44 , 
		       45  , 46 , 47, 48 , 49 , 50 , 51 ,  53 , 54 , 55};
  
};


void CRTDataQuality::ResetVar()
{

  run             = -9999;
  subrun          = -9999;
  event           = -9999;
  date            = -9999; 
  nFlashes        = -9999; 
  nMatchedFlashes = -9999; 
  nTracks         = -9999; 
  nMatchedTracks  = -9999; 
  
  for (size_t i= 0; i < 73; i++ )
    {
      nCRThits[i]  = 0;
      AvgPe[i]     = 0.;
    }
   reaoutTime = -9999.;
   CRTFlashMatchTimeDiff.clear();
}

CRTDataQuality::CRTDataQuality(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    data_labelhit_(p.get<std::string>("data_labelhit_","crthitcorr")),
    data_label_DAQHeader_(p.get<std::string>("data_label_DAQHeader_","daq")),  
    data_label_TPCtrack_(p.get<std::string>("data_label_TPCtrack_","pandora")),  
    data_label_flash_(p.get<std::string>("data_label_flash_","simpleFlashCosmic")),  
    data_label_match_(p.get<std::string>("data_label_match_","crttrackmatch")),  
    minT_(p.get<double>("minT_",-1500000.)),
    maxT_(p.get<double>("maxT_", 3500000.)),
    flash_match_timecut(p.get<double>("flash_match_timecut",3.0)),  // in us
    min_length_TPCtrack(p.get<double>("min_length_TPCtrack",20.0)),  // in cm
    fTimeZeroOffset(p.get<int>("fTimeZeroOffset",60000)),
    applyTimeOffSet_(p.get<bool>("applyTimeOffSet",true)),
    verbose_(p.get<bool>("verbose",true))
    // More initializers here.    
{
}

void CRTDataQuality::analyze(art::Event const & evt)
{
  ResetVar();


  // This will be useful for the date
  art::Timestamp        evtTime  = evt.time();
  long int timeInNsSec = (long int) evtTime.value () ;
  date = (int) (timeInNsSec >>  32 ); // Convert ns time stamp into second time stamp

  
  // -------------------------------------- Check if all we need is valid -----------------------------------------
  //Get GPS Time
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
  double evt_timeGPS_nsec = (double)evtTimeGPS.timeLow();  


  //get CRTHits
  art::Handle< std::vector<crt::CRTHit> > rawHandle_hit;
  evt.getByLabel(data_labelhit_, rawHandle_hit); 
  //check to make sure the data we asked for is valid
  if(!rawHandle_hit.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTHits " << " in module " << data_labelhit_ << std::endl;
    std::cout << std::endl;
    return;
  } // This should throw an exception, but ok WRONG!


  // get TPC Track List
  //check whether tracks exist
  bool iTrackAssnCRT = false; 
  art::Handle< std::vector<recob::Track>  > trackListHandle; 
  std::vector<art::Ptr<recob::Track> >  tracklist;
  if (evt.getByLabel(data_label_TPCtrack_,trackListHandle))
    {
      if(trackListHandle.isValid())
	{ 
	  iTrackAssnCRT = true;
	  art::fill_ptr_vector(tracklist, trackListHandle);
	}else{
	std::cout << "no data product found for crt-track matching" << std::endl;
      }
    }else {
    std::cout << "no data product found for crt-track matching" << std::endl;
  }
    
  //check whether tzeros exist
  bool iT0crt = false;
  art::Handle< std::vector<anab::T0> > rawHandle_Tzero;
  evt.getByLabel(data_label_match_, rawHandle_Tzero);
  if(rawHandle_Tzero.isValid()) {
    // grab T0 objects associated with tracks    
    iT0crt=true;
    //    art::FindMany<anab::T0> trk_t0C_assn_v(trackListHandle, evt, data_label_match_);
  }
  else {
    std::cout << "no data product found for crt-track matching" << std::endl;
  }


  
  //get Optical Flash Collection
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlash;
  evt.getByLabel(data_label_flash_, rawHandle_OpFlash);  
  std::vector<recob::OpFlash> const& OpFlashCollection(*rawHandle_OpFlash);
  if(verbose_){ 
    std::cout<<"  OpFlashCollection.size()  "<<OpFlashCollection.size()<<std::endl; 
  }  //get Optical Flash

  // -------------------------------------- End of Check if all we need is valid -----------------------------------------

  run    = evt.run() ;
  subrun = evt.subRun() ;
  event  = evt.event() ;
  // Build FEB to ArrayIndex Conversion
  std::map<int, int> febIndexConversion;
  for (size_t iM = 0; iM<73 ; iM++ ) febIndexConversion[febIndex[iM] ] = iM;
 
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle_hit);
  
  // Loop over CRT Hits, keep them if they're within some readout time
  for(size_t j = 0; j < CRTHitCollection.size(); j++) {   
    crt::CRTHit my_CRTHit = CRTHitCollection[j];
    // Calculate the reaout time (time we're considering the CRT hit in)
    reaoutTime = maxT_ - minT_;
    // Calculate the time of this CRT Hit
    double thisHitTime =(double)my_CRTHit.ts0_ns;
    double offset      =  ((double)fTimeZeroOffset - evt_timeGPS_nsec);
    if ( applyTimeOffSet_ ) thisHitTime += offset; 

    hHitTime->Fill(thisHitTime);
    // If the hit time is not within the readout time, skip
    if (thisHitTime < minT_ ) continue;
    if (thisHitTime > maxT_ ) continue;
    hHitTimeAfterCut->Fill(thisHitTime);
    // Every hit is composed by 2 FEBs: the horizontal and vertical one
    auto febID_v = my_CRTHit.feb_id;
    for (auto const febID : febID_v) {
      int key   = (int) febID;
      int index = febIndexConversion[key];
     
      nCRThits[ index ]++ ;
      AvgPe   [ index ] += my_CRTHit.peshit ;
      hModule_X[ index ]->Fill(my_CRTHit.x_pos);
      hModule_Y[ index ]->Fill(my_CRTHit.y_pos);
      hModule_Z[ index ]->Fill(my_CRTHit.z_pos);
    }         
        
  }//Loop on CRT hits
  

  for (size_t i = 0; i < 73; i++)
    {
      if (nCRThits[i]){  AvgPe[ i ] /= (float)nCRThits[i];}
      else AvgPe[ i ] = -999.;
    } 

  CRTFlashMatchTimeDiff.clear();
  calculateFlashMatch(CRTHitCollection, OpFlashCollection, evt_timeGPS_nsec,nFlashes, nMatchedFlashes, CRTFlashMatchTimeDiff);

  if (iT0crt&&iTrackAssnCRT)  
    {
      art::FindMany<anab::T0> trk_t0C_assn_v(trackListHandle, evt, data_label_match_);
      calculateTrackMatch(trackListHandle, trk_t0C_assn_v, tracklist.size(), nTracks ,nMatchedTracks );
    }
  // load the tree for this event
  fTree->Fill();
 

}



void  CRTDataQuality::calculateFlashMatch( std::vector<crt::CRTHit> const CRTHitCollection,    
					   std::vector<recob::OpFlash> const OpFlashCollection,
					   double event_timeGPS_ns, int &nFlashes, int &nMatchedFlashes, std::vector<double> &CRTFlashMatchTimeDiff)
{
  
  int ntotFlash = OpFlashCollection.size();
  if (ntotFlash>100) ntotFlash=100;
  nMatchedFlashes=0;
  nFlashes=ntotFlash;
  int iFlashM[100] = {0};
  
  
  // Loop over CRT Hits, keep them if they're within some readout time
  for(size_t j = 0; j < CRTHitCollection.size(); j++) {   
    crt::CRTHit my_CRTHit = CRTHitCollection[j];
    // Calculate the time of this CRT Hit
    
    double thisHitTime =(double)my_CRTHit.ts0_ns;
    double offset      =  ((double)fTimeZeroOffset - event_timeGPS_ns);
    if ( applyTimeOffSet_ ) thisHitTime += offset; 
    // Loop over PMT flashes
    float bestdiff =  flash_match_timecut;
    //    thisHitTime += (double)fTimeZeroOffset;
    double signedbestdiff;
    int ibestflash = -1;
    for(int i = 0; i != ntotFlash; i++) {
      if (iFlashM[i]==0) {
	recob::OpFlash my_OpFlash = OpFlashCollection[i];
	auto Timeflash = my_OpFlash.Time(); //in us from trigger time
	float thisdiff = fabs(0.001*thisHitTime- Timeflash);
	if (thisdiff<bestdiff) {
	  bestdiff=thisdiff;
	  signedbestdiff = Timeflash-0.001*thisHitTime;
	  ibestflash=i;
	}
      }
    }// end loop on PMT flashes    
    if (ibestflash>=0) {
      iFlashM[ibestflash]=1;
      hFlashDiff->Fill(signedbestdiff);
      nMatchedFlashes++;
      CRTFlashMatchTimeDiff.push_back(signedbestdiff);
    }       
  }//Loop on CRT hits
  if(verbose_)std::cout << ntotFlash << " " << nMatchedFlashes << std::endl;
  
}

void  CRTDataQuality::calculateTrackMatch(  art::Handle< std::vector<recob::Track>  > trackListHandle, 
					     art::FindMany<anab::T0> trk_t0C_assn_v,
					    int Tcollsize, int &nTracks ,int &nMatchedTracks )
{
    
      
  nMatchedTracks=0;
  int nTPCtracks = Tcollsize;
  nTracks=0;
  for(int j = 0; j < nTPCtracks; j++) {
    art::Ptr<recob::Track> ptrack(trackListHandle, j);
    const recob::Track& track = *ptrack;
    
    if (track.Length()>min_length_TPCtrack) {
      nTracks++;
      const std::vector<const anab::T0*>& T0_v = trk_t0C_assn_v.at(j);
      if (T0_v.size()==1) { 
	auto t0 = T0_v.at(0);
	//	  int tzerotime =t0->Time();	  
	int plane =t0->TriggerBits();  // this variable is not filled until MCC9.1, defaults to 0
	double dca = t0->TriggerConfidence();
	if (plane==3) hTrackDCA3->Fill(dca);
	else  hTrackDCA->Fill(dca);
	nMatchedTracks++;
      }  // if we have t0tags available
    }  // if length is at least 20 cm
  } // end loop over tracks
  //  std::cout << nTPCtracks << " " << denom << " " << numer << std::endl;

}



void CRTDataQuality::beginJob()
{
  // Implementation of optional member function here.
  
  fTree = tfs->make<TTree>("CRTDataQuality","analysis tree");
  fTree->Branch("run"       ,&run       ,"run/I"   );
  fTree->Branch("subrun"    ,&subrun    ,"subrun/I");
  fTree->Branch("event"     ,&event     ,"event/I" );
  fTree->Branch("date"      ,&date      ,"date/I"  );
  fTree->Branch("febIndex"  ,febIndex   ,"febIndex[73]/I");
  fTree->Branch("AvgPe"     ,AvgPe      ,"AvgPe[73]/D"   );
  fTree->Branch("nFlashes"         ,&nFlashes        ,"nFlashes/I"         );
  fTree->Branch("nMatchedFlashes"  ,&nMatchedFlashes ,"nMatchedFlashes/I"  );
  fTree->Branch("nTracks"          ,&nTracks         ,"nTracks/I"          );
  fTree->Branch("nMatchedTracks"   ,&nMatchedTracks  ,"nMatchedTracks/I"   );
  fTree->Branch("nCRThits"  ,nCRThits   ,"nCRThits[73]/I");
  fTree->Branch("febIndex"  ,febIndex   ,"febIndex[73]/I");
  fTree->Branch("AvgPe"     ,AvgPe      ,"AvgPe[73]/D"   );
  fTree->Branch("CRTFlashMatchTimeDiff","std::vector<double>",&CRTFlashMatchTimeDiff);

  /*
  fTree->Branch("AvgCountX" ,AvgCountX  ,"AvgCountX[73]/D");
  fTree->Branch("AvgCountY" ,AvgCountY  ,"AvgCountY[73]/D");
  fTree->Branch("AvgCountZ" ,AvgCountZ  ,"AvgCountZ[73]/D");
  fTree->Branch("MaxCountX" ,MaxCountX  ,"MaxCountX[73]/D");
  fTree->Branch("MaxCountY" ,MaxCountY  ,"MaxCountY[73]/D");
  fTree->Branch("MaxCountZ" ,MaxCountZ  ,"MaxCountZ[73]/D");
  fTree->Branch("MaxPositX" ,MaxPositX  ,"MaxPositX[73]/D");
  fTree->Branch("MaxPositY" ,MaxPositY  ,"MaxPositY[73]/D");
  fTree->Branch("MaxPositZ" ,MaxPositZ  ,"MaxPositZ[73]/D");
  */
  fTree->Branch("reaoutTime",&reaoutTime,"reaoutTime/D"  );
  
  hHitTime         = tfs1->make<TH1D>("hHitTime"        ,"CRT Hit Time; time [ns]; ",2000, -10000000,10000000);
  hHitTimeAfterCut = tfs1->make<TH1D>("hHitTimeAfterCut","CRT Hit Time; time [ns]; ",2000, -10000000,10000000);
  hFlashDiff = tfs1->make<TH1D>("hFlashDiff"," ",1000, -5.,5.);
  hFlashDiff->GetXaxis()->SetTitle("Flash time - CRT Hit time (us)");
  hTrackDCA = tfs1->make<TH1D>("hTrackDCA"," ",200,0.,50.);
  hTrackDCA->GetXaxis()->SetTitle("DCA of TPC track and CRT hit side and bottom planes (cm))");
  hTrackDCA3 = tfs1->make<TH1D>("hTrackDCA3"," ",200,0.,50.);
  hTrackDCA3->GetXaxis()->SetTitle("DCA of TPC track and CRT hit top plane only (cm))");
  
  for (size_t i = 0; i<73; i++)
    {
      std::string feb = std::to_string(febIndex[i]);
      std::string x = "hModule_X" + feb;
      std::string y = "hModule_Y" + feb;
      std::string z = "hModule_Z" + feb;
      hModule_X[i]  = tfs1->make<TH1D>(x.c_str() ,"CRT Hit X Position; X [cm]; ",2000, -1000, 1000);
      hModule_Y[i]  = tfs1->make<TH1D>(y.c_str() ,"CRT Hit Y Position; Y [cm]; ",2000, -1000, 1000); 
      hModule_Z[i]  = tfs1->make<TH1D>(z.c_str() ,"CRT Hit Z Position; Z [cm]; ",2000, -500 , 1500);

    }

}
 
void CRTDataQuality::endJob()
{  

}


DEFINE_ART_MODULE(CRTDataQuality)


