///////////////////////////////////////////////////////////////////////
// Class:       T0recoCRT
// Module Type: producer
// File:        T0recoCRT_module.cc
// Description: It generates and associates a T0 object to TPC Tracks ussing CRT
// Generated at Fri Jun 29 09:34:08 2018 by David Lorca Galindo using artmod
// from cetpkgsupport v1_14_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
//#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// data-products                                                                         
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
  class T0recoCRT;
}

class crt::T0recoCRT : public art::EDProducer {
public:
  explicit T0recoCRT(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  T0recoCRT(T0recoCRT const &) = delete;
  T0recoCRT(T0recoCRT &&) = delete;
  T0recoCRT & operator = (T0recoCRT const &) = delete;
  T0recoCRT & operator = (T0recoCRT &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  //from Chris                                                                                                                                            
  void SortTrackPoints (const recob::Track& track, std::vector<recob::Track::Point_t>& sorted_trk);


private:

  // Declare member data here.
  art::ServiceHandle<art::TFileService> tfs;

  std::string  data_labelCRTtrack_;
  std::string  data_labelCRThit_;
  std::string  data_label_flash_;
  std::string  data_label_DAQHeader_;
  std::string  data_label_TPCTrack_;
  int fHardDelay_;
  int fCRTT0off_;
  unsigned int fShowerCut_;
  double fThetaCut_;
  double fPhiCut_;
  int fGPSMatchW_;
  double fvdrift_;
  int storeAsn_;
  int verbose_;


  
};


crt::T0recoCRT::T0recoCRT(fhicl::ParameterSet const & p)
 :
// Initialize member data here.
  data_labelCRTtrack_(p.get<std::string>("data_labelCRTtrack")),
  data_labelCRThit_(p.get<std::string>("data_labelCRThit")),
  data_label_flash_(p.get<std::string>("data_label_flash_")),
  data_label_DAQHeader_(p.get<std::string>("data_label_DAQHeader_")),
  data_label_TPCTrack_(p.get<std::string>("data_label_TPCTrack_")),
  fHardDelay_(p.get<int>("fHardDelay",40000)),
  fCRTT0off_(p.get<int>("fCRTT0off",69000)),
  fShowerCut_(p.get<int>("fShowerCut",40)),
  fThetaCut_(p.get<double>("fThetaCut",5.12)),
  fPhiCut_(p.get<double>("fPhiCut",12.12)),
  fGPSMatchW_(p.get<int>("fGPSMatchW",500)),
  fvdrift_(p.get<double>("fvdrift",0.111436)),
  storeAsn_(p.get<int>("storeAsn",1)),
  verbose_(p.get<int>("verbose"))  // ,
{
  // Call appropriate produces<>() functions here.

  produces< std::vector< anab::T0 > >();
  produces<art::Assns<recob::Track, recob::OpFlash> >();
  produces< art::Assns <recob::Track, anab::T0> >();
  produces< art::Assns <recob::Track, crt::CRTTrack > >();
  
}

void crt::T0recoCRT::produce(art::Event & evt)
{
  // Implementation of required member function here.

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
  evt.getByLabel(data_labelCRThit_, rawHandle_CRThit); //                                                                                                       //check to make sure the data we asked for is valid                                                                                                
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
    getchar();
  }                                                                                                                                        
  //get CRTTracks 


  // produce data-product to be filled within module
  std::unique_ptr< std::vector<anab::T0> > T0_v(new std::vector<anab::T0>);
  std::unique_ptr< art::Assns <recob::Track, anab::T0> >       trk_t0_assn_v   ( new art::Assns<recob::Track, anab::T0>);
  std::unique_ptr< art::Assns<recob::Track, recob::OpFlash> > trk_flash_assn_v (new art::Assns<recob::Track, recob::OpFlash>);
  std::unique_ptr< art::Assns <recob::Track, crt::CRTTrack > > trk_crttrack_assn_v( new art::Assns<recob::Track, crt::CRTTrack > );

  art::PtrMaker<recob::Track> trackPtrMaker(evt, rawHandle_TPCtrack.id());
  art::PtrMaker<recob::OpFlash> flashPtrMaker(evt, rawHandle_OpFlash.id());
  art::PtrMaker<crt::CRTTrack> crttrackPtrMaker(evt, rawHandle_CRTtrack.id());
  //art::PtrMaker<anab::T0> t0PtrMaker(evt, *this);  
  art::PtrMaker<anab::T0> t0PtrMaker(evt);  
  
  if( CRTTrackCollection.size()<fShowerCut_) {//A0 cut in showers
    
  for(std::vector<int>::size_type i = 0; i != TPCTrackCollection.size(); i++) {//A 
    
    recob::Track my_TPCTrack = TPCTrackCollection[i];

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
      
      if( (ThetaDiffABS<fThetaCut_) && (PhiDiffABS<fPhiCut_) ){//C //Angular Cut
	
	int CRTTrack_T1_nsec = my_CRTTrack.ts1_ns + fHardDelay_; 
	int CRTTrack_T0_nsec = my_CRTTrack.ts0_ns + fCRTT0off_;
	
	for(std::vector<int>::size_type j = 0; j != OpFlashCollection.size(); j++) {//D //look for flash in time with CRTTrack
	  
	  recob::OpFlash my_OpFlash = OpFlashCollection[j];
	  auto Timeflash = my_OpFlash.Time(); //in us from trigger time                 
	  auto Timeflash_ns = (Timeflash * 1000);
	  auto Timeflash_ns_GPS = evt_timeGPS_nsec + (Timeflash * 1000);
	  
	  int TdiffT1_nsec = Timeflash_ns - CRTTrack_T1_nsec;	
  	  int TdiffT0_nsec = Timeflash_ns_GPS - CRTTrack_T0_nsec;
	  int TdiffT0_nsecABS = std::abs(TdiffT0_nsec);
	  
	  if(TdiffT0_nsecABS<fGPSMatchW_){//E : cut in GPS Match <1us
	  
	    //for this flash, make T0 object and associations
	    double crtT = CRTTrack_T0_nsec;
	    double dT = CRTTrack_T0_nsec - Timeflash_ns_GPS;	   
	    
	    anab::T0 t0(crtT, 0, 1, 1, dT);
	    T0_v->emplace_back(t0);
	    
	    //make pointers and associations
	    art::Ptr<recob::Track> trackptr = trackPtrMaker(i);
	    art::Ptr<recob::OpFlash> flashptr = flashPtrMaker(j);
	    art::Ptr<crt::CRTTrack> crttrackptr = crttrackPtrMaker(k);
	    art::Ptr<anab::T0> t0ptr = t0PtrMaker(T0_v->size()-1);
	    
	    trk_flash_assn_v->addSingle(trackptr,flashptr);
	    trk_crttrack_assn_v->addSingle(trackptr,crttrackptr);
	    trk_t0_assn_v->addSingle(trackptr,t0ptr);
	    
	    
	    if(verbose_!=0 ){
	      std::cout.precision(19);  
	      std::cout<<"CRT & TPC Tracks match in angular cuts"<<std::endl;
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
	      getchar();    
	    }
	    
	  }//E
	  
	}//D
	
      }//C
      
    }//B
    
  }//A
  
  }//A0 cut in showers 
  
  if(storeAsn_ == 1){
    evt.put(std::move(T0_v));
    evt.put(std::move(trk_t0_assn_v));
    evt.put(std::move(trk_flash_assn_v));
    evt.put(std::move(trk_crttrack_assn_v));
  }
  
}

void crt::T0recoCRT::beginJob()
{
  // Implementation of optional member function here.
}

void crt::T0recoCRT::endJob()
{
  // Implementation of optional member function here.
}

void crt::T0recoCRT::SortTrackPoints(const recob::Track& track, std::vector<recob::Track::Point_t>& sorted_trk)
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


DEFINE_ART_MODULE(crt::T0recoCRT)
