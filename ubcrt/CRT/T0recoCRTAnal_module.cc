////////////////////////////////////////////////////////////////////////
// Class:       T0recoCRTAnal
// Module Type: analyzer
// File:        T0recoCRTAnal_module.cc
//
// Generated at Tue Jul 10 08:05:45 2018 by David Lorca Galindo using artmod
// from cetpkgsupport v1_14_01.
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

//data-products
#include "lardataobj/RecoBase/Track.h"
//#include "Track.h"
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
  class T0recoCRTAnal;
}

class crt::T0recoCRTAnal : public art::EDAnalyzer {
public:
  explicit T0recoCRTAnal(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  T0recoCRTAnal(T0recoCRTAnal const &) = delete;
  T0recoCRTAnal(T0recoCRTAnal &&) = delete;
  T0recoCRTAnal & operator = (T0recoCRTAnal const &) = delete;
  T0recoCRTAnal & operator = (T0recoCRTAnal &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void SortTrackPoints (const recob::Track& track, std::vector<recob::Track::Point_t>& sorted_trk);

private:

  // Declare member data here.
  art::ServiceHandle<art::TFileService> tfs;

  std::string  data_labelCRTtrack_;
  std::string  data_labelCRThit_;
  std::string  data_label_flash_;
  std::string  data_label_DAQHeader_;
  std::string  data_label_TPCTrack_;
  std::string  data_label_T0reco_;
  int fHardDelay_;
  int fCRTT0off_;
  double fvdrift_;
  int storeAsn_;
  int verbose_;

  TH1F* hDiffT_CRT_T0;
  TH1F* hDiffT_CRT_T0Flash;
  TH1F* hDiffT_T0Flash;
  TH1F* hGeoMatch;
  TH1F* hFlashCount;
  TH1F* hFlashCountClone;

};


crt::T0recoCRTAnal::T0recoCRTAnal(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  data_labelCRTtrack_(p.get<std::string>("data_labelCRTtrack")),
  data_labelCRThit_(p.get<std::string>("data_labelCRThit")),
  data_label_flash_(p.get<std::string>("data_label_flash_")),
  data_label_DAQHeader_(p.get<std::string>("data_label_DAQHeader_")),
  data_label_TPCTrack_(p.get<std::string>("data_label_TPCTrack_")),
  data_label_T0reco_(p.get<std::string>("data_label_T0reco_")),
  fHardDelay_(p.get<int>("fHardDelay",40000)),
  fCRTT0off_(p.get<int>("fCRTT0off",69000)),
  fvdrift_(p.get<double>("fvdrift",0.111436)),
  storeAsn_(p.get<int>("storeAsn",1)),
  verbose_(p.get<int>("verbose"))  // ,
 // More initializers here.
{}

void crt::T0recoCRTAnal::analyze(art::Event const & evt)
{



  //get DAQHeader for GPS time                                                                                                                                 
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
  auto evt_timeGPS_sec = evtTimeGPS.timeHigh();
  auto evt_timeGPS_nsec = evtTimeGPS.timeLow();
  if(verbose_!=0){
    std::cout<<"Evt:GPS_Time: sec, nsec :"<<evt_timeGPS_sec <<" "<< evt_timeGPS_nsec <<std::endl;
  }
  //get DAQHeader for GPS time                                                                                                                                  
  
  
  //get TPC Tracks  
  art::Handle< std::vector<recob::Track> > rawHandle_TPCtrack;
  evt.getByLabel(data_label_TPCTrack_, rawHandle_TPCtrack);
  //check to make sure the data we asked for is valid                                                                                                       
  
  if(!rawHandle_TPCtrack.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " recob::Track " << " in module " << data_label_TPCTrack_ << std::endl;
    std::cout << std::endl;
    return;
  }
  //get better access to the data	    
  
  std::vector<recob::Track> const& TPCTrackCollection(*rawHandle_TPCtrack);
  if(verbose_!=0){
    std::cout<<"  TPCTrackCollection.size()  "<<TPCTrackCollection.size()<<std::endl;
    // getchar();			   
    
  }
  //get TPCTracks                                                                                                                                               
  
  //get Optical Flash			   
  
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlash;
  evt.getByLabel(data_label_flash_, rawHandle_OpFlash);
  std::vector<recob::OpFlash> const& OpFlashCollection(*rawHandle_OpFlash);
  if(verbose_!=0){
    std::cout<<"  OpFlashCollection.size()  "<<OpFlashCollection.size()<<std::endl;
  }
  //get Optical Flash 
  
  //get CRTHits				   
  art::Handle< std::vector<crt::CRTHit> > rawHandle_CRThit;
  evt.getByLabel(data_labelCRThit_, rawHandle_CRThit); //                                                                                                       
  //check to make sure the data we asked for is valid	
  
  if(!rawHandle_CRThit.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTHits " << " in module " << data_labelCRThit_ << std::endl;
    std::cout << std::endl;
    return;
  }
  
  //get better access to the data      
  
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle_CRThit);
  if(verbose_!=0){
    std::cout<<"  CRTHitCollection.size()  "<<CRTHitCollection.size()<<std::endl;
    //  getchar();			  
  }
  //get CRTHits                                                                                                                                                 
  
  //get CRTTracks			    
  
  art::Handle< std::vector<crt::CRTTrack> > rawHandle_CRTtrack;
  evt.getByLabel(data_labelCRTtrack_, rawHandle_CRTtrack);
  //check to make sure the data we asked for is valid	
  
  if(!rawHandle_CRTtrack.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTTracks " << " in module " << data_labelCRTtrack_ << std::endl;
    std::cout << std::endl;
    return;
  }
  //get better access to the data                                                                                                                             
  
  std::vector<crt::CRTTrack> const& CRTTrackCollection(*rawHandle_CRTtrack);
  if(verbose_!=0){
    std::cout<<"  CRTTrackCollection.size()  "<<CRTTrackCollection.size()<<std::endl;
    // getchar();			      
    
  }
  //get CRTTracks                                                                                                                                               
  
  // grab T0 objects associated with tracks    
  art::FindMany<anab::T0> trk_t0_assn_v(rawHandle_TPCtrack, evt, data_label_T0reco_); //objeto, evento, label
  
  // grab flashes associated with tracks 
  art::FindMany<recob::OpFlash> trk_flash_assn_v(rawHandle_TPCtrack, evt, data_label_T0reco_ );
  
  if( CRTTrackCollection.size()<40) {//A0 cut in showers

    for(std::vector<int>::size_type i = 0; i != TPCTrackCollection.size(); i++) {//A                                  
      
      recob::Track my_TPCTrack = TPCTrackCollection[i];
      //auto TPCTrackLength =  my_TPCTrack.Length();
      
      std::vector<recob::Track::Point_t> sorted_trk;
      SortTrackPoints(my_TPCTrack,sorted_trk);
      
      auto const& top    = sorted_trk.at(0);
      auto const& bottom = sorted_trk.at(sorted_trk.size() - 1);
      
      double TT[] = {top.X(),top.Y(),top.Z()};
      double BT[] = {bottom.X(),bottom.Y(),bottom.Z()};
      double VT[] = {BT[0] - TT[0], BT[1] - TT[1], BT[2] - TT[2]};
      
      double TPCTheta = crt::auxfunctions::CalTheta(VT[0],VT[1],VT[2]);
      double TPCPhi = crt::auxfunctions::CalPhi(VT[0],VT[1],VT[2]);
      
      if(verbose_!=0){
	std::cout<<"TPC Track Director Vector: ("<<VT[0]<<","<<VT[1]<<","<<VT[2]<<")"<<std::endl;
	std::cout<<"TPC Theta: "<<TPCTheta<<std::endl;
	std::cout<<"TPC Phi: "<<TPCPhi<<std::endl;
      }
      
      for(std::vector<int>::size_type k = 0; k != CRTTrackCollection.size(); k++) {//B                                                                           
	
	crt::CRTTrack my_CRTTrack = CRTTrackCollection[k];
	
	double TC[3]={0,0,0}, BC[3]={0,0,0}, VC[3]={0,0,0};
	
	if(my_CRTTrack.z1_pos>my_CRTTrack.z2_pos){
	  TC[0] = my_CRTTrack.x1_pos;
	  TC[1] = my_CRTTrack.y1_pos;
	  TC[2] = my_CRTTrack.z1_pos;
	  BC[0] = my_CRTTrack.x2_pos;
	  BC[1] = my_CRTTrack.y2_pos;
	  BC[2] = my_CRTTrack.z2_pos;
	}
	
	if(my_CRTTrack.z2_pos>my_CRTTrack.z1_pos){
	  TC[0] = my_CRTTrack.x2_pos;
	  TC[1] = my_CRTTrack.y2_pos;
	  TC[2] = my_CRTTrack.z2_pos;
	  BC[0] = my_CRTTrack.x1_pos;
	  BC[1] = my_CRTTrack.y1_pos;
	  BC[2] = my_CRTTrack.z1_pos;
	}
	
	VC[0] = BC[0] - TC[0];
	VC[1] = BC[1] - TC[1];
	VC[2] = BC[2] - TC[2];
	
	double CRTTheta = crt::auxfunctions::CalTheta(VC[0],VC[1],VC[2]);
	double CRTPhi = crt::auxfunctions::CalPhi(VC[0],VC[1],VC[2]);
	
	if(verbose_!=0){
	  std::cout<<"CRT Track Director Vector: ("<<VC[0]<<","<<VC[1]<<","<<VC[2]<<")"<<std::endl;
	  std::cout<<"CRT Theta: "<<CRTTheta<<std::endl;
	  std::cout<<"CRT Phi: "<<CRTTheta<<std::endl;
	  getchar();
	  
	}
	
	double ThetaDiff = TPCTheta-CRTTheta;
	double ThetaDiffABS = fabs(ThetaDiff);
	double PhiDiff = TPCPhi-CRTPhi;
	double PhiDiffABS = fabs(PhiDiff);
	
	if( (ThetaDiffABS<5.1) && (PhiDiffABS<12.1) ){//C //Angular Cut                                                                                          
	  
	  int CRTTrack_T1_nsec = my_CRTTrack.ts1_ns + fHardDelay_;
	  int CRTTrack_T0_nsec = my_CRTTrack.ts0_ns + fCRTT0off_;	  

	  int countFlash = 0;
	  for(std::vector<int>::size_type j = 0; j != OpFlashCollection.size(); j++) {//D //look for flash in time with CRTTrack                             	  
	    recob::OpFlash my_OpFlash = OpFlashCollection[j];
	    auto Timeflash = my_OpFlash.Time(); //in us from trigger time                                                                                   
	    auto Timeflash_ns = (Timeflash * 1000);
	    auto Timeflash_ns_GPS = evt_timeGPS_nsec + (Timeflash * 1000);
	    
	    int TdiffT1_nsec = Timeflash_ns - CRTTrack_T1_nsec;
	    int TdiffT0_nsec = Timeflash_ns_GPS - CRTTrack_T0_nsec;	    
	    int TdiffT0_nsecABS = std::abs(TdiffT0_nsec);


	    /* if( ((TdiffT1_nsec)>310)  &&  ((TdiffT1_nsec)<570) ) countFlash++;
	    if( ((TdiffT1_nsec)> -8700)  &&  ((TdiffT1_nsec)< -8620) ) countFlash++;
	    if( ((TdiffT1_nsec)>8400)  &&  ((TdiffT1_nsec)<8480) ) countFlash++;
	    if( ((TdiffT1_nsec)>10500)  &&  ((TdiffT1_nsec)<10580) ) countFlash++;
	    
	    if(countFlash==1)hGeoMatch->Fill(TdiffT1_nsec);
	    */
  
	    //if( ((TdiffT1_nsec)>300)  &&  ((TdiffT1_nsec)<600) ){//E:: cut in time T1  
	    if(TdiffT0_nsecABS<250){//E : cut in GPS Match   	    

	      //make associations between TPC Track, CRT Track, Flash and T0.


	    if(verbose_!=0){
	      std::cout.precision(19);
	      std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
	      std::cout<<"CRT & TPC Tracks match in angular and time"<<std::endl;
	      std::cout<<"TPC Track: "<<i<<" out of "<<CRTTrackCollection.size()<<std::endl;
	      std::cout<<"     "<<std::endl;
	      std::cout<<"Matched Flash["<<j<<"]: "<<Timeflash_ns<<" ns w.r.t. trigger"<<std::endl;
	      std::cout<<"CRTTrack Time(T1corr["<<k<<"]: "<<CRTTrack_T1_nsec<<" ns w.r.t. trigger"<<std::endl;
	      std::cout<<"CRTTrack Time T1["<<k<<"]: "<<my_CRTTrack.ts1_ns<<" ns w.r.t. trigger"<<std::endl;
	      std::cout<<"hardDelay: "<<fHardDelay_<<std::endl;
	      std::cout<<"     "<<std::endl;
	      std::cout<<"Time Difference: "<<TdiffT1_nsec<<" ns"<<std::endl;
	      std::cout<<"     "<<std::endl;
	      std::cout<<"TPCTrack_Theta: "<<TPCTheta<<"    TPCTrack_Phi: "<<TPCPhi<<std::endl;
	      std::cout<<"CRTTrack_Theta: "<<CRTTheta<<"    CRTTrack_Phi: "<<CRTPhi<<std::endl;
	      std::cout<<"ThetaDiffABS: "<<ThetaDiffABS<<std::endl;
	      std::cout<<"PhiDiffABS: "<<PhiDiffABS<<std::endl;
	      //getchar();
	    }
	    
            
            const std::vector<const anab::T0*>& T0_v = trk_t0_assn_v.at(i);
	    
            const std::vector<const recob::OpFlash*>& flash_v = trk_flash_assn_v.at(i);
	    
	    if( (flash_v.size() != 0) && (T0_v.size() != 0)  ){//F                                                                                             
	      
              auto t0 = T0_v.at(0);
              auto TimeT0 = t0->Time();
              auto TimeT0_ns = (TimeT0 * 1000);
	      
              auto T0flash = flash_v.at(0);
              auto T0Timeflash = T0flash->Time(); //in us from trigger time
              auto T0Timeflash_ns = (T0Timeflash * 1000);

              int T0diff = TimeT0_ns - CRTTrack_T1_nsec;
              int T0Flashdiff = T0Timeflash_ns - CRTTrack_T1_nsec;

	      hDiffT_CRT_T0->Fill(T0diff);
	      hDiffT_CRT_T0Flash->Fill(T0Flashdiff);
	      hDiffT_T0Flash->Fill(TimeT0_ns - T0Timeflash_ns);
	      
	      if(verbose_!=0){
		std::cout<<"     "<<std::endl;
		std::cout<<"AND for this ONE I also have a T0 and a Flash"<<std::endl;
		std::cout<<"     "<<std::endl;
		std::cout<<"T0 time is: "          <<TimeT0_ns<<" ns"<<std::endl;
		std::cout<<"T0 associated flash: " <<T0Timeflash_ns<<" ns"<<std::endl;
		
		
		std::cout<<"Tdiff T0-CRTtrack:      " <<T0diff<<" ns"<<std::endl;
		std::cout<<"Tdiff T0Flash-CRTtrack: " <<T0Flashdiff<<" ns"<<std::endl;
		std::cout<<"     "<<std::endl;
		std::cout<<"AND for this ONE I also have a T0 and a Flash"<<std::endl;
		
		for(std::vector<int>::size_type a = 0; a != OpFlashCollection.size(); a++){
		  recob::OpFlash my_OpFlashA = OpFlashCollection[a];
		  std::cout<<"Flash["<<a<<"] Time: "<<my_OpFlashA.Time() * 1000<<"ns w.r.t. trigger "<<std::endl;	      
		}
		std::cout<<"     "<<std::endl;
		std::cout<<"     "<<std::endl;
		std::cout<<"     "<<std::endl;
		getchar();
	      }
	      
            }//F
	    else{//F2
	      
	      if(verbose_!=0){
		std::cout<<"For this match I do not have a T0 object"<<std::endl;
		getchar();
	      }

	    }//F2                                                                                                                                               
	    
	    }//E  //cut in beam T1                                                                                                                                               
            
	  }//D                                                                                                                                                   
	  hFlashCount->Fill(countFlash);
	  hFlashCountClone->Fill(countFlash);
	  
	}//C                                                                                                                                                     
	
      }//B                                                                                                                                                       

    }//A 

  }//A0 cut in showers



}

void crt::T0recoCRTAnal::beginJob()
{

  hDiffT_CRT_T0 = tfs->make<TH1F>("hDiffT_CRT_T0","hDiffT_CRT_T0",5000,-20000,20000);
  hDiffT_CRT_T0->GetXaxis()->SetTitle("T0Time - CRTTrack_Time (ns)");
  hDiffT_CRT_T0->GetYaxis()->SetTitle("Entries/bin");

  hDiffT_CRT_T0Flash = tfs->make<TH1F>("hDiffT_CRT_T0Flash","hDiffT_CRT_T0Flash",5000,-20000,20000);
  hDiffT_CRT_T0Flash->GetXaxis()->SetTitle("T0Time_Flash - CRTTrack_Time (ns)");
  hDiffT_CRT_T0Flash->GetYaxis()->SetTitle("Entries/bin");

  hDiffT_T0Flash = tfs->make<TH1F>("hDiffT_T0Flash","hDiffT_T0Flash",5000,-20000,20000);
  hDiffT_T0Flash->GetXaxis()->SetTitle("T0Time - T0Time_Flash (ns)");
  hDiffT_T0Flash->GetYaxis()->SetTitle("Entries/bin");

  hGeoMatch = tfs->make<TH1F>("hGeoMatch","hGeoMatch",5000,-20000,20000);
  hGeoMatch->GetXaxis()->SetTitle("FlashTime - CRTTrack_Time (ns)");
  hGeoMatch->GetYaxis()->SetTitle("Entries/bin");

  hFlashCount = tfs->make<TH1F>("hFlashCount","hFlashCount",6,-1,5);
  hFlashCount->GetXaxis()->SetTitle("Flashes satisfying GeoCut");
  hFlashCount->GetYaxis()->SetTitle("Entries/bin");

  hFlashCountClone = tfs->make<TH1F>("hFlashCountFrac","hFlashCountFrac",6,-1,5);
  hFlashCountClone->GetXaxis()->SetTitle("Flashes satisfying GeoCut");
  hFlashCountClone->GetYaxis()->SetTitle("Fraction of events");


}

void crt::T0recoCRTAnal::endJob()
{
  // Implementation of optional member function here.

  Double_t norm = hFlashCountClone->GetEntries();
  hFlashCountClone->Scale(1/norm);


}

void crt::T0recoCRTAnal::SortTrackPoints(const recob::Track& track, std::vector<recob::Track::Point_t>& sorted_trk)
{

  sorted_trk.clear();

  auto const&N = track.NumberTrajectoryPoints();
  auto const&start = track.LocationAtPoint(0);
  auto const&end   = track.LocationAtPoint( N - 1 );

  if (start.Y() > end.Y()){
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( track.LocationAtPoint(i) );
  }

  else {
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( track.LocationAtPoint( N - i - 1) );
  }
}

DEFINE_ART_MODULE(crt::T0recoCRTAnal)
