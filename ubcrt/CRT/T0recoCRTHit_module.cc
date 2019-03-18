//////////////////////////////////////////////////////////////////////
// Class:       T0recoCRTHit
// Module Type: producer
// File:        T0recoCRTHit_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include <artdaq-core/Data/Fragment.hh>

#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTTzero.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3S.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "Math/SMatrix.h"
#include "Math/Functions.h"
#include "TDatime.h"
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <typeinfo>

class T0recoCRTHit : public art::EDProducer {
public:
  explicit T0recoCRTHit(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  T0recoCRTHit(T0recoCRTHit const &) = delete;
  T0recoCRTHit(T0recoCRTHit &&) = delete;
  T0recoCRTHit & operator = (T0recoCRTHit const &) = delete;
  T0recoCRTHit & operator = (T0recoCRTHit &&) = delete;

 // Required functions.
  void produce(art::Event & e);
  
  // Selected optional functions.
  void beginJob();
  void endJob();

private:

  art::ServiceHandle<art::TFileService> tfs;
  
  uint32_t fEvtNum; //Number of current event                       
  uint32_t frunNum;                //Run Number taken from event  
  uint32_t fsubRunNum;             //Subrun Number taken from event         
  std::string  data_label_TPCtrack_;
  std::string  data_label_CRTtzero_;
  std::string  data_label_CRThit_;
  std::string  data_label_flash_;
  std::string  data_label_DAQHeader_;
  int fHardDelay;
  int fTimeZeroOffset;
  int fTimeSelect;
  int fMatchCutTop;
  int fMatchCut;
  float fDriftVel;
  bool fstoreAssn;
  bool fverbose;
  bool fIsMC;
  //alignment params
  float fAlignBotX;
  float fAlignBotY;
  float fAlignBotZ;
  float fAlignAnodeX;
  float fAlignAnodeY;
  float fAlignAnodeZ;
  float fAlignCathX;
  float fAlignCathY;
  float fAlignCathZ;
  float fAlignTopX;
  float fAlignTopY;
  float fAlignTopZ;

  //histograms

  TH1F* hTrackLength;
  TH1F* hOpeningAngle;

  // TH1F* hDistTwo;
  TH1F* hDistOne;
  TH1F* hDistAn;
  TH1F* hDistCa;
  TH1F* hDistNone;
  TH1F* hTimeBest;

  /*
  TH1F* hTrackLengthACPT;
  TH1F* hOpeningAngleACPT;
  TH1F* hThetaAll;
  TH1F* hThetaACPT;
  TH1F* hThetaGood;
  TH1F* hThetaGoodA;
  TH1F* hThetaGoodY;
  TH1F* hThetaGoodN;
  TH1F* hThetaAnodeY;
  TH1F* hThetaAnodeN;
  TH1F* hThetaCathY;
  TH1F* hThetaCathN;
  TH1F* hPhiAll;
  TH1F* hPhiACPT;
  TH1F* hPhiGood;
  TH1F* hPhiGoodA;
  TH1F* hPhiGoodY;
  TH1F* hPhiGoodN;
  TH1F* hPhiAnodeY;
  TH1F* hPhiAnodeN;
  TH1F* hPhiCathY;
  TH1F* hPhiCathN;
  TH2F* hTvsPACPT;
  TH2F* hTvsPGood;
  TH2F* hTvsPGoodA;
  TH2F* hTvsPGoodY;
  TH2F* hTvsPGoodN;
  TH2F* hTvsPAnodeY;
  TH2F* hTvsPAnodeN;
  TH2F* hTvsPCathY;
  TH2F* hTvsPCathN;
  TH1F* hDistTrue[4];
  TH1F* hDistWrong[4];
  TH1F* hDistTrueAll;
  TH1F* hDistWrongAll;
  TH1F* hDistNone;
  TH1F* hDistTrueA;
  TH1F* hDistWrongA;
  TH1F* hDistNoneA;
  TH1F* hPlaneClosest;
  TH1F* hPlaneCorr;
  TH1F* hPlaneWrong;
  TH1F* hDistMistake;
  TH1F* hTimeDiffCorr; 
  TH1F* hTimeDiffWrong;
  TH1F* hTimeDiffNone;
  */

  
};


T0recoCRTHit::T0recoCRTHit(fhicl::ParameterSet const & p)
  : 
    data_label_TPCtrack_(p.get<std::string>("data_label_TPCtrack")),
    data_label_CRTtzero_(p.get<std::string>("data_label_CRTtzero")),
    data_label_CRThit_(p.get<std::string>("data_label_CRThit")),
    data_label_flash_(p.get<std::string>("data_label_flash")),
    data_label_DAQHeader_(p.get<std::string>("data_label_DAQHeader")),
    fHardDelay(p.get<int>("HardDelay",40000)),
    fTimeZeroOffset(p.get<int>("TimeZeroOffset",60000)),
    fTimeSelect(p.get<int>("TimeSelect",0)),
    fMatchCutTop(p.get<int>("MatchCutTop",50)),
    fMatchCut(p.get<int>("MatchCut",25)),
    fDriftVel(p.get<float>("DriftVel",0.11436)),   // cm/us
    fstoreAssn(p.get<bool>("storeAssn",true)),
    fverbose(p.get<bool>("verbose",false)),
    fIsMC(p.get<bool>("IsMC",false)),
    fAlignBotX(p.get<float>("AlignBotX",0.0)),
    fAlignBotY(p.get<float>("AlignBotY",0.0)),
    fAlignBotZ(p.get<float>("AlignBotZ",0.0)),
    fAlignAnodeX(p.get<float>("AlignAnodeX",0.0)),
    fAlignAnodeY(p.get<float>("AlignAnodeY",0.0)),
    fAlignAnodeZ(p.get<float>("AlignAnodeZ",0.0)),
    fAlignCathX(p.get<float>("AlignCathX",0.0)),
    fAlignCathY(p.get<float>("AlignCathY",0.0)),
    fAlignCathZ(p.get<float>("AlignCathZ",0.0)),
    fAlignTopX(p.get<float>("AlignTopX",0.0)),
    fAlignTopY(p.get<float>("AlignTopY",0.0)),
    fAlignTopZ(p.get<float>("AlignTopZ",0.0))
{

  // Call appropriate produces<>() functions here.
  produces< std::vector< anab::T0 > >();
  //  produces<art::Assns<recob::Track, recob::OpFlash> >();
  produces< art::Assns <recob::Track, anab::T0> >();
  produces< art::Assns <recob::Track, crt::CRTTzero > >();

}


void T0recoCRTHit::produce(art::Event & evt)
{
  // Implementation of required member function here.
  
  frunNum    = evt.run();
  fsubRunNum = evt.subRun();
  fEvtNum = evt.event();
  
  // get TPC Track List 
  art::Handle< std::vector<recob::Track>  > trackListHandle;
  std::vector<art::Ptr<recob::Track> >  tracklist;

  if (evt.getByLabel(data_label_TPCtrack_,trackListHandle))
      art::fill_ptr_vector(tracklist, trackListHandle);

  // get any t0 information for these tracks from ACPT module
  art::FindMany<anab::T0> trk_t0_assn_v(trackListHandle, evt, "t0reco" );
  //  art::FindMany<recob::OpFlash> trk_flash_assn_v(trackListHandle, evt, "t0reco" );

  
  std::unique_ptr< std::vector<anab::T0> > T0_v(new std::vector<anab::T0>);
  std::unique_ptr< art::Assns <recob::Track, anab::T0> >       trk_t0_assn_v_new   ( new art::Assns<recob::Track, anab::T0>);
  //  std::unique_ptr< art::Assns<recob::Track, recob::OpFlash> > trk_flash_assn_v (new art::Assns<recob::Track, recob::OpFlash>);
  std::unique_ptr< art::Assns <recob::Track, crt::CRTTzero > > trk_crttzero_assn_v( new art::Assns<recob::Track, crt::CRTTzero > );


  double evt_timeGPS_nsec=0.0;   double evt_timeGPS_sec=0.0;
  if (!fIsMC) {  // this is data

    //check to make sure the data we asked for is valid 
    //get DAQ Header                                                                  
    art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;  
    evt.getByLabel(data_label_DAQHeader_, rawHandle_DAQHeader);
    
    if(!rawHandle_DAQHeader.isValid()){
      std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
		<< ", event " << evt.event() << " has zero"
		<< " DAQHeaderTimeUBooNE  " << " in with label " << data_label_DAQHeader_ << std::endl;    
      evt.put(std::move(T0_v));
      evt.put(std::move(trk_t0_assn_v_new));
      evt.put(std::move(trk_crttzero_assn_v));
      return;
    }
    
    raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
    art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();  
    evt_timeGPS_sec = evtTimeGPS.timeHigh();
    evt_timeGPS_nsec = (double)evtTimeGPS.timeLow();
    // art::Timestamp evtTimeNTP = my_DAQHeader.ntp_time();
    // double evt_timeNTP_sec = evtTimeNTP.timeHigh();
    // double evt_timeNTP_nsec = (double)evtTimeNTP.timeLow();
    // double timstp_diff = std::abs(evt_timeGPS_nsec - evt_timeNTP_nsec);
  }// end if data

  if(fverbose){
    std::cout<< "Run:  "<<frunNum << "   subRun: " <<fsubRunNum<<std::endl;
    std::cout<<"event: "<<fEvtNum <<std::endl;
    std::cout.precision(19);
    std::cout<<"  GPS time second:  "<<evt_timeGPS_sec<<std::endl;
    std::cout<<"  GPS time nano_second:  "<<evt_timeGPS_nsec<<std::endl;
    // std::cout<<"  NTP time second:  "<<evt_timeNTP_sec<<std::endl;    
    // std::cout<<"  NTP time nano_second:  "<<evt_timeNTP_nsec<<std::endl;
    // std::cout<<"  event time second:  "<<evt_time_sec<<std::endl;
    // std::cout<<"  event time nano_second:  "<<evt_time_nsec<<std::endl;
    // std::cout<<"  difference between GPS and NTP:  "<<evt_timeGPS_nsec - evt_timeNTP_nsec<<" ns"<<std::endl;
    // std::cout<<"  ABS difference between GPS and NTP:  "<<timstp_diff<<" ns"<<std::endl;
    // if( (evt_time_sec==evt_timeGPS_sec) && (evt_time_nsec==evt_timeGPS_nsec))  std::cout<<" Event time type is: GPS  "<<std::endl;
    // if( (evt_time_sec==evt_timeNTP_sec) && (evt_time_nsec==evt_timeNTP_nsec))  std::cout<<" Event time type is: NTP  "<<std::endl;
  }// end if verbose  
  
   

  // fetch tzeros 
  art::Handle< std::vector<crt::CRTTzero> > rawHandletzero;
  evt.getByLabel(data_label_CRTtzero_, rawHandletzero); //what is the product instance name? no BernZMQ
  //check to make sure the data we asked for is valid                                           
  if(!rawHandletzero.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTTzeros " << " in module " << data_label_CRTtzero_ << std::endl;
    std::cout << std::endl;
    evt.put(std::move(T0_v));
    evt.put(std::move(trk_t0_assn_v_new));
    evt.put(std::move(trk_crttzero_assn_v));
    return;
  }
  std::vector<art::Ptr<crt::CRTTzero> > tzerolist;
  if (evt.getByLabel(data_label_CRTtzero_,rawHandletzero))
    art::fill_ptr_vector(tzerolist, rawHandletzero);
  art::FindManyP<crt::CRTHit> fmht(rawHandletzero, evt, data_label_CRTtzero_);
  

  // produce data-product to be filled within module
  art::PtrMaker<recob::Track> trackPtrMaker(evt, trackListHandle.id());
  //art::PtrMaker<recob::OpFlash> flashPtrMaker(evt, rawHandle_OpFlash.id());
  art::PtrMaker<crt::CRTTzero> crttzPtrMaker(evt, rawHandletzero.id());
  //  art::PtrMaker<anab::T0> t0PtrMaker(evt, *this);  
  art::PtrMaker<anab::T0> t0PtrMaker(evt);  

  //get Optical Flash
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlash;
  evt.getByLabel(data_label_flash_, rawHandle_OpFlash);  
  std::vector<recob::OpFlash> const& OpFlashCollection(*rawHandle_OpFlash);
  if(fverbose){ 
    std::cout<<"  OpFlashCollection.size()  "<<OpFlashCollection.size()<<std::endl; 
  }  //get Optical Flash
  

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
	
	// if (fTimeSelect==0) // for EXT BNB data
	//     xshift = ((double)tzerolist[tzIter]->ts0_ns - (double)evt_timeGPS_nsec)*vdrift);
	// else //fTimeSelect_==1 for BNB data
	float diff;
	if (fTimeSelect==1) diff= fabs(0.001*(tzerolist[tzIter]->ts1_ns+fHardDelay)-Timeflash);
	else diff = fabs(0.001*(tzerolist[tzIter]->ts0_ns+fTimeZeroOffset-(double)evt_timeGPS_nsec)-Timeflash);
	if (diff<min_deltat) { min_deltat=diff; best_time_match=tzIter;}
      } // loop over tzeros     
      if (best_time_match>=0) {
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
   
  
// Set up space charge map
  //Spacecharge services provider 
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
   
  double const vdrift = fDriftVel*0.001;  // in cm/ns
  // const detinfo::DetectorProperties *_detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // double const vdrift =  _detprop->DriftVelocity();  
  // std::cout << "drift velocity is " << vdrift << std::endl;

  // loop over tracks  
  if (tracklist.size()>0) {
    for (size_t trkIter = 0; trkIter < tracklist.size(); ++trkIter) {
      
      //    int TrackGood=0;
      float dist_besthit = 99999999.; int plane_besthit=111;
      double time_besthit=999999999.; 
      double theta,phi,trklen;
      // fetch track length, start and end points, theta and phi
      auto startP=tracklist[trkIter]->Start();      
      auto endP=tracklist[trkIter]->End();      
      trklen = tracklist[trkIter]->Length();
      theta=tracklist[trkIter]->Theta();      phi=tracklist[trkIter]->Phi();
      if (phi>0) {theta=3.14159-theta; phi=phi-3.14159;}
      
      // get track directional cosines
      auto trackCosStart = tracklist[trkIter]->StartDirection();
      auto trackCosEnd = tracklist[trkIter]->EndDirection();

      double opang = trackCosStart.X()*trackCosEnd.X() +  trackCosStart.Y()*trackCosEnd.Y() + 
	trackCosStart.Z()*trackCosEnd.Z();
      
      //reject tracks that are too short and bend too much
      if (trklen>20 && opang>0.95)  {
	
	if (fverbose) {
	  std::cout << "Event " << evt.event() <<  " Track " << trkIter << " cos(opening angle) " << 
	    opang << "  track length " << trklen << std::endl;
	  std::cout << "        theta " << theta << " phi " << phi << std::endl;
	}
	
	int besttz=-1; 
	//loop over CRT tzeros
	if (tzerolist.size()>0) {   
	  for(size_t tzIter = 0; tzIter < tzerolist.size(); ++tzIter){   
	    
	    //calculate track shift in x for the time  of this CRT hit
	    double xshift = 0.0;
	    if (fTimeSelect==0) // for EXT BNB data
	      xshift = ((double)tzerolist[tzIter]->ts0_ns +fTimeZeroOffset- (double)evt_timeGPS_nsec)*vdrift;
	    else //fTimeSelect==1 for BNB data
	      xshift = ((double)(tzerolist[tzIter]->ts1_ns) + fHardDelay)*vdrift;	      
	    // verify that track is within the fiducial volume after this time shift
	    double test1 = startP.X()-xshift; double test2 = endP.X()-xshift;
	    if (test1>-10. && test1<270. && test2>-10 && test2<270.) {
	      
	    //  loop over CRT hits for this tzero
	      std::vector<art::Ptr<crt::CRTHit> > hitlist=fmht.at(tzIter);
	      if (hitlist.size()>0) {
		for (size_t ah = 0; ah< hitlist.size(); ++ah){	
		  
		  //adjust CRT hit positions for rough alignment correction 
		  //   temporary until the CRT geometry is updated with the precise positions from the laser survey
		  double crt_x=hitlist[ah]->x_pos;
		  double crt_y=hitlist[ah]->y_pos;
		  double crt_z=hitlist[ah]->z_pos;
		  
		  //		  int crt_plane = (hitlist[ah]->plane)%10;
		  /*
		  if (crt_plane==0) {
		    crt_x+=fAlignBotX;
		    crt_y+=fAlignBotY;
		    crt_z+=fAlignBotZ;
		  }
		  else if (crt_plane==1) {
		    crt_x+=fAlignAnodeX;
		    crt_y+=fAlignAnodeY;
		    crt_z+=fAlignAnodeZ;
		  }
		  else if (crt_plane==2) {
		    crt_x+=fAlignCathX;
		    crt_y+=fAlignCathY;
		    crt_z+=fAlignCathZ;
		  }
		  else if (crt_plane==3) {
		    crt_x+=fAlignTopX;
		    crt_y+=fAlignTopY;
		    crt_z+=fAlignTopZ;
		  }
		  */
		  TVector3 CRTpoint(crt_x,crt_y,crt_z);
		  
		  //calculate track shift in x for the time  of this CRT hit
		  if (fTimeSelect==0) // for EXT BNB data
		    xshift = ((double)hitlist[ah]->ts0_ns + fTimeZeroOffset - evt_timeGPS_nsec)*vdrift;
		  else //fTimeSelect_==1 for BNB data
		    xshift = ((double)(hitlist[ah]->ts1_ns) + fHardDelay)*vdrift;	      
		  
		  //  Correct start and end point for space charge with this tzero
		  geo::Point_t newStartP = startP; geo::Point_t newEndP = endP;
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
		    
		    // std::cout << fPosOffsets.X() << " " <<   fPosOffsets.Y() << " " <<  fPosOffsets.Z() << std::endl;
		  }
		  else {
		    newStartP.SetX(startP.X()-xshift); newEndP.SetX(endP.X()-xshift);
		  }

		  TVector3 trackstart(newStartP.X(),newStartP.Y(),newStartP.Z());
		  TVector3 trackend(newEndP.X(),newEndP.Y(),newEndP.Z());
		  
		  // calculate the distance of closest approach (DCA) of track to CRT hit
		  TVector3 denom = trackend-trackstart;
		  TVector3 a =CRTpoint-trackstart;  TVector3 b=CRTpoint-trackend;
		  TVector3 numer = a.Cross(b);
		  double dca = numer.Mag()/denom.Mag();
		  
		  if (dca<dist_besthit) {
		    besttz = tzIter;
		    dist_besthit=dca;
		    plane_besthit=(hitlist[ah]->plane)%10;
		    //	    art::Ptr<recob::OpFlash> flashptr = flashPtrMaker(j);
		    //	    art::Ptr<crt::CRTHit> crthitptr = crthitPtrMaker(k);
		    
		    if (fTimeSelect==0)
		      time_besthit = (double)hitlist[ah]->ts0_ns + fTimeZeroOffset- evt_timeGPS_nsec;
		    else  time_besthit=hitlist[ah]->ts1_ns+ fHardDelay; 
		  }
		} // loop over CRT hits
	      }// if hits 
	    }  // track fid vol check
	  } //loop over tzeros
	}  //if tzeros
	
	// fill histograms and creates associations	
	if (dist_besthit<200) {
	  //	  std::vector<art::Ptr<crt::CRTHit> > hitlist=fmht.at(besttz);
	  // int icountplanes =1;
	  // if (hitlist.size()>0) {
	  //   for (size_t ah = 0; ah< hitlist.size(); ++ah){	
	  //     if (hitlist[ah]->plane!=plane_besthit) icountplanes=2;
	  //   }
	  // }
	  hDistOne->Fill(dist_besthit);

		  
	  //calculate track shift in x for the time  of this CRT hit
	  double xshift = time_besthit*vdrift;
	  TVector3 trackstart(startP.X()-xshift,startP.Y(),startP.Z());
	  TVector3 trackend(endP.X()-xshift,endP.Y(),endP.Z());
	  
	  if (plane_besthit==1) {
	    if ((trackstart.X()>-10 && trackstart.X()<2)  || (trackend.X()>-10 && trackend.X()<2) )
	      hDistAn->Fill(dist_besthit); 
	  }
	  else if (plane_besthit==2) {
	    if ((trackstart.X()>255 && trackstart.X()<268)  || (trackend.X()>255 && trackend.X()<268) )
	      hDistCa->Fill(dist_besthit); 
	  }	  
	  if ((dist_besthit<fMatchCut) || (dist_besthit<fMatchCutTop && plane_besthit==3)) {
	    //	    double dT =0.0;
	    // args are (time, triggertype, triggerbits, ?, trigger confidence)
	    anab::T0 thist0(0.001*time_besthit, 2, 1, 1, dist_besthit);
	    T0_v->emplace_back(thist0);
	    
	    //make pointers and associations
	    art::Ptr<recob::Track> trackptr = trackPtrMaker(trkIter);
	    //	    art::Ptr<recob::OpFlash> flashptr = flashPtrMaker(j);
	    art::Ptr<crt::CRTTzero> crttzeroptr = crttzPtrMaker(besttz);
	    art::Ptr<anab::T0> t0ptr = t0PtrMaker(T0_v->size()-1);
	    
	    //	    trk_flash_assn_v->addSingle(trackptr,flashptr);
	    trk_crttzero_assn_v->addSingle(trackptr,crttzeroptr);
	    trk_t0_assn_v_new->addSingle(trackptr,t0ptr);
	  }
	  else hDistNone->Fill(100.);
	  //
	  hTimeBest->Fill(time_besthit*0.001);
	} // if hit is close	  
      } // if good track
      
      
      //    if (TrackGood) {       hThetaGood->Fill(theta);      hPhiGood->Fill(phi); hTvsPGood->Fill(phi,theta);}
      hTrackLength->Fill(trklen);
      hOpeningAngle->Fill(opang);
      //    hThetaAll->Fill(theta);
      //hPhiAll->Fill(phi);    
      
    } //loop over tracks
  }  // if tracks
  
  if(fstoreAssn){
    evt.put(std::move(T0_v));
    evt.put(std::move(trk_t0_assn_v_new));
    //    evt.put(std::move(trk_flash_assn_v));
    evt.put(std::move(trk_crttzero_assn_v));
  }
  

}

void T0recoCRTHit::beginJob()
{
  // Implementation of optional member function here.

  // Create histograms 
  hTrackLength = tfs->make<TH1F>("hTrackLength","hTrackLength",175,0.,350.);
  hTrackLength->GetXaxis()->SetTitle("Track Length (cm)");
  //  
  hOpeningAngle = tfs->make<TH1F>("hOpeningAngle","hOpeningAngle",50,0.,1.);
  hOpeningAngle->GetXaxis()->SetTitle("cosine opening angle");
  //
  
  hDistOne = tfs->make<TH1F>("hDistOne","hDistOne",200,0.,200.);
  hDistOne->GetXaxis()->SetTitle("distance (cm)");
  hDistAn = tfs->make<TH1F>("hDistAn","hDistAn",200,0.,200.);
  hDistAn->GetXaxis()->SetTitle("distance (cm)");
  hDistCa = tfs->make<TH1F>("hDistCa","hDistCa",200,0.,200.);
  hDistCa->GetXaxis()->SetTitle("distance (cm)");
  // hDistTwo= tfs->make<TH1F>("hDistTwo","hDistTwo",200,0.,200.);    
  // hDistTwo->GetXaxis()->SetTitle("distance (cm)");
  hDistNone= tfs->make<TH1F>("hDistNone","hDistNone",200,0.,200.);    
  hDistNone->GetXaxis()->SetTitle("distance (cm)");
  hTimeBest = tfs->make<TH1F>("hTimeBest","hTimeBest",200,0.,200.);
  hTimeBest->GetXaxis()->SetTitle("time (us)");

  /*
  hThetaAll = tfs->make<TH1F>("hThetaAll","hThetaAll",100,0.,3.1416);
  hThetaAll->GetXaxis()->SetTitle("theta (radians)");
  hThetaACPT= tfs->make<TH1F>("hThetaACPT","hThetaACPT",100,0.,3.1416);
  hThetaACPT->GetXaxis()->SetTitle("theta (radians)");
  hThetaGood = tfs->make<TH1F>("hThetaGood","hThetaGood",100,0.,3.1416);
  hThetaGood->GetXaxis()->SetTitle("theta (radians)");
  hThetaGoodA = tfs->make<TH1F>("hThetaGoodA","hThetaGoodA",100,0.,3.1416);
  hThetaGoodA->GetXaxis()->SetTitle("theta (radians)");
  hThetaGoodY = tfs->make<TH1F>("hThetaGoodY","hThetaGoodY",100,0.,3.1416);
  hThetaGoodY->GetXaxis()->SetTitle("theta (radians)");
  hThetaGoodN = tfs->make<TH1F>("hThetaGoodN","hThetaGoodN",100,0.,3.1416);
  hThetaGoodN->GetXaxis()->SetTitle("theta (radians)");
  hThetaAnodeY = tfs->make<TH1F>("hThetaAnodeY","hThetaAnodeY",100,0.,3.1416);
  hThetaAnodeY->GetXaxis()->SetTitle("theta (radians)");
  hThetaAnodeN = tfs->make<TH1F>("hThetaAnodeN","hThetaAnodeN",100,0.,3.1416);
  hThetaAnodeN->GetXaxis()->SetTitle("theta (radians)");
  hThetaCathY = tfs->make<TH1F>("hThetaCathY","hThetaCathY",100,0.,3.1416);
  hThetaCathY->GetXaxis()->SetTitle("theta (radians)");
  hThetaCathN = tfs->make<TH1F>("hThetaCathN","hThetaCathN",100,0.,3.1416);
  hThetaCathN->GetXaxis()->SetTitle("theta (radians)");
  hPhiAll= tfs->make<TH1F>("hPhiAll","hPhiAll",100,-3.1416,0.);
  hPhiAll->GetXaxis()->SetTitle("phi (radians)");
  hPhiACPT = tfs->make<TH1F>("hPhiACPT","hPhiACPT",100,-3.1416,0.);
  hPhiACPT->GetXaxis()->SetTitle("phi (radians)");
  hPhiGood= tfs->make<TH1F>("hPhiGood","hPhiGood",100,-3.1416,0.);
  hPhiGood->GetXaxis()->SetTitle("phi (radians)");
  hPhiGoodA= tfs->make<TH1F>("hPhiGoodA","hPhiGoodA",100,-3.1416,0.);
  hPhiGoodA->GetXaxis()->SetTitle("phi (radians)");
  hPhiGoodY= tfs->make<TH1F>("hPhiGoodY","hPhiGoodY",100,-3.1416,0.);
  hPhiGoodY->GetXaxis()->SetTitle("phi (radians)");
  hPhiGoodN= tfs->make<TH1F>("hPhiGoodN","hPhiGoodN",100,-3.1416,0.);
  hPhiGoodN->GetXaxis()->SetTitle("phi (radians)");
  hPhiAnodeY= tfs->make<TH1F>("hPhiAnodeY","hPhiAnodeY",100,-3.1416,0.);
  hPhiAnodeY->GetXaxis()->SetTitle("phi (radians)");
  hPhiAnodeN= tfs->make<TH1F>("hPhiAnodeN","hPhiAnodeN",100,-3.1416,0.);
  hPhiAnodeN->GetXaxis()->SetTitle("phi (radians)");
  hPhiCathY= tfs->make<TH1F>("hPhiCathY","hPhiCathY",100,-3.1416,0.);
  hPhiCathY->GetXaxis()->SetTitle("phi (radians)");
  hPhiCathN= tfs->make<TH1F>("hPhiCathN","hPhiCathN",100,-3.1416,0.);
  hPhiCathN->GetXaxis()->SetTitle("phi (radians)");
  hTvsPACPT= tfs->make<TH2F>("hTvsPACPT","hTvsPACPT",50,-3.1416,0.,50,0.,3.1416);
  hTvsPACPT->GetXaxis()->SetTitle("phi (radians)");
  hTvsPACPT->GetYaxis()->SetTitle("theta (radians)");
  hTvsPGood= tfs->make<TH2F>("hTvsPGood","hTvsPGood",50,-3.1416,0.,50,0.,3.1416);
  hTvsPGood->GetXaxis()->SetTitle("phi (radians)");
  hTvsPGood->GetYaxis()->SetTitle("theta (radians)");
  hTvsPGoodA= tfs->make<TH2F>("hTvsPGoodA","hTvsPGoodA",50,-3.1416,0.,50,0.,3.1416);
  hTvsPGoodA->GetXaxis()->SetTitle("phi (radians)");
  hTvsPGoodA->GetYaxis()->SetTitle("theta (radians)");
  hTvsPGoodY= tfs->make<TH2F>("hTvsPGoodY","hTvsPGoodY",50,-3.1416,0.,50,0.,3.1416);
  hTvsPGoodY->GetXaxis()->SetTitle("phi (radians)");
  hTvsPGoodY->GetYaxis()->SetTitle("theta (radians)");
  hTvsPGoodN= tfs->make<TH2F>("hTvsPGoodN","hTvsPGoodN",50,-3.1416,0.,50,0.,3.1416);
  hTvsPGoodN->GetXaxis()->SetTitle("phi (radians)");
  hTvsPGoodN->GetYaxis()->SetTitle("theta (radians)");
  hTvsPAnodeY= tfs->make<TH2F>("hTvsPAnodeY","hTvsPAnodeY",50,-3.1416,0.,50,0.,3.1416);
  hTvsPAnodeY->GetXaxis()->SetTitle("phi (radians)");
  hTvsPAnodeY->GetYaxis()->SetTitle("theta (radians)");
  hTvsPAnodeN= tfs->make<TH2F>("hTvsPAnodeN","hTvsPAnodeN",50,-3.1416,0.,50,0.,3.1416);
  hTvsPAnodeN->GetXaxis()->SetTitle("phi (radians)");
  hTvsPAnodeN->GetYaxis()->SetTitle("theta (radians)");
  hTvsPCathY= tfs->make<TH2F>("hTvsPCathY","hTvsPCathY",50,-3.1416,0.,50,0.,3.1416);
  hTvsPCathY->GetXaxis()->SetTitle("phi (radians)");
  hTvsPCathY->GetYaxis()->SetTitle("theta (radians)");
  hTvsPCathN= tfs->make<TH2F>("hTvsPCathN","hTvsPCathN",50,-3.1416,0.,50,0.,3.1416);
  hTvsPCathN->GetXaxis()->SetTitle("phi (radians)");
  hTvsPCathN->GetYaxis()->SetTitle("theta (radians)");


  //
  for (int i=0;i<4;++i) {
    TString label1 = Form("hDistTrue%1d",i);
    hDistTrue[i] = tfs->make<TH1F>(label1,label1,100,0.,100.);
    TString label2 = Form("hDistWrong%1d",i);
    hDistWrong[i]= tfs->make<TH1F>(label2,label2,100,0.,100.);    
  }
  //
  hDistTrueAll = tfs->make<TH1F>("hDistTrueAll","hDistTrueAll",100,0.,100.);
  hDistWrongAll= tfs->make<TH1F>("hDistWrongAll","hDistWrongAll",100,0.,100.);    
  hDistNone= tfs->make<TH1F>("hDistNone","hDistNone",100,0.,100.);    
  //
  hDistTrueA = tfs->make<TH1F>("hDistTrueA","hDistTrueA",100,0.,100.);
  hDistWrongA =  tfs->make<TH1F>("hDistWrongA","hDistWrongA",100,0.,100.);    
  hDistNoneA = tfs->make<TH1F>("hDistNoneA","hDistNoneA",100,0.,100.);    
  //
  hPlaneClosest = tfs->make<TH1F>("PlaneClosest","PlaneClosest",4,-0.5,3.5);
  hPlaneCorr = tfs->make<TH1F>("PlaneCorr","PlaneCorr",4,-0.5,3.5);
  hPlaneWrong = tfs->make<TH1F>("PlaneWrong","PlaneWrong",4,-0.5,3.5);
  //
  hDistMistake = tfs->make<TH1F>("hDistMistake","hDistMistake",100,0.,10000.);
  //
  hTimeDiffCorr = tfs->make<TH1F>("hTimeDiffCorr","hTimeDiffCorr",500,-25.,25.);
  hTimeDiffWrong = tfs->make<TH1F>("hTimeDiffWrong","hTimeDiffWrong",500,-25.,25.);
  hTimeDiffNone = tfs->make<TH1F>("hTimeDiffNone","hTimeDiffNone",500,-25.,25.);
  */

    /*  
  hFlashTimeDis = tfs->make<TH1F>("hFlashTimDis","hFlashTimDis",2000,-5,25);
  hFlashTimeDis->GetXaxis()->SetTitle("Flash Time w.r.t. trigger (us)");
  hFlashTimeDis->GetYaxis()->SetTitle("Entries/bin");  

  hFlashTimeDis_b0 = tfs->make<TH1F>("hFlashTimDis_b0","hFlashTimDis_b0",2000,-5,25);
  hFlashTimeDis_b0->GetXaxis()->SetTitle("Flash_bo Time w.r.t. trigger (us)");
  hFlashTimeDis_b0->GetYaxis()->SetTitle("Entries/bin");

  hTFvsTH_t1 = tfs->make<TH1F>("hBeamMatching","hBeamMatching",500,-1000,1000);
  hTFvsTH_t1->GetXaxis()->SetTitle("Flash Time w.r.t. trigger - CRTHit Time_t1 (ns)");
  hTFvsTH_t1->GetYaxis()->SetTitle("Entries/bin");

  hTFvsTH_t1_2d = tfs->make<TH2F>("hBeamMatching2","hBeamMatching2",6,-3,3,500,0,1000);
  hTFvsTH_t1_2d->GetXaxis()->SetTitle("Flash Time - CRTHit Time (s)");
  hTFvsTH_t1_2d->GetYaxis()->SetTitle("Flash Time w.r.t Trigger - CRTHit_Time_t1 (ns)");
  hTFvsTH_t1_2d->SetOption("COLZ"); 

  hTFvsTH_t0 = tfs->make<TH1F>("hGPSMatching","GPSMatching",2000,0,2000000);
  hTFvsTH_t0->GetXaxis()->SetTitle("Flash_Time_GPS - CRTHit_Time_T0 (ns)");
  hTFvsTH_t0->GetYaxis()->SetTitle("Entries/bin");

  hTFvsTH_t0_2d = tfs->make<TH2F>("hGPSMatching2","hGPSMatching2",6,-3,3,2000,0,2000000);
  hTFvsTH_t0_2d->GetXaxis()->SetTitle("Flash Time - CRTHit Time (s)");
  hTFvsTH_t0_2d->GetYaxis()->SetTitle("Flasf_Time_GPS - CRTHit_Time_T0 (ns)");
  hTFvsTH_t0_2d->SetOption("COLZ"); 

  hTFvsTH_t0_t1 = tfs->make<TH2F>("hGPSBeamMatching","hGPSBeamMatching",2000,0,2000000,500,0,1000);
  hTFvsTH_t0_t1->GetXaxis()->SetTitle("Flash_Time_GPS - CRTHit_Time_T0 (ns)");
  hTFvsTH_t0_t1->GetYaxis()->SetTitle("Flash Time w.r.t Trigger - CRTHit_Time_t1 (ns)");
  hTFvsTH_t0_t1->SetOption("COLZ"); 

  hTFvsTH_plane_t0 = tfs->make<TH2F>("hGPSBeamMatchingPlane","hGPSBeamMatchingPlane",4,0,4, 4000,0,2000000);
  hTFvsTH_plane_t0->GetXaxis()->SetTitle("CRT plane (0=bottom, 1=FT, 2=Pipe, 3=Top))");
  hTFvsTH_plane_t0->GetYaxis()->SetTitle("Flash Time_GPS - CRTHit Time_t0 (ns)");
  hTFvsTH_plane_t0->SetOption("COLZ"); 


  hNFlavsNHit = tfs->make<TH2F>("hNFlavsNHit","hNFlavsNHit",30,0,30,100,0,300);
  hNFlavsNHit->GetXaxis()->SetTitle("Number of Flashes per event");
  hNFlavsNHit->GetYaxis()->SetTitle("Number of CRT Hits per event");
  hNFlavsNHit->GetZaxis()->SetTitle("Entries/bin");
  hNFlavsNHit->SetOption("COLZ");

  hNHitperFla = tfs->make<TH1F>("hNHitperFla","hNHitperFla",205,-5,200);
  hNHitperFla->GetXaxis()->SetTitle("N^{o} of CRTHits per Flash");
  hNHitperFla->GetYaxis()->SetTitle("Entries/bin");

  hNHitperFla0 = tfs->make<TH1F>("hNHitperFlaBot","hNHitperFlaBot",205,-5,200);
  hNHitperFla0->GetXaxis()->SetTitle("N^{o} of CRTHits per Flash in Bottom");
  hNHitperFla0->GetYaxis()->SetTitle("Entries/bin");

  hNHitperFla1 = tfs->make<TH1F>("hNHitperFlaFT","hNHitperFlaFT",205,-5,200);
  hNHitperFla1->GetXaxis()->SetTitle("N^{o} of CRTHits per Flash in FT");
  hNHitperFla1->GetYaxis()->SetTitle("Entries/bin");

  hNHitperFla2 = tfs->make<TH1F>("hNHitperFlaPipe","hNHitperFlaPipe",205,-5,200);
  hNHitperFla2->GetXaxis()->SetTitle("N^{o} of CRTHits per Flash in Pipe");
  hNHitperFla2->GetYaxis()->SetTitle("Entries/bin");

  hNHitperFla3 = tfs->make<TH1F>("hNHitperFlaTop","hNHitperFlaTop",205,-5,200);
  hNHitperFla3->GetXaxis()->SetTitle("N^{o} of CRTHits per Flash in Top");
  hNHitperFla3->GetYaxis()->SetTitle("Entries/bin");
	
  hNHitperFla2D = tfs->make<TH2F>("hNHitperEvtPlane","hNHitperEvtPlane",4,0,4,205,-5,200);
  hNHitperFla2D->GetXaxis()->SetTitle("CRT plane (0=bottom, 1=FT, 2=Pipe, 3=Top))");
  hNHitperFla2D->GetYaxis()->SetTitle("N^{o} of CRTHits in event");
  hNHitperFla2D->SetOption("COLZ"); 


  hNTraperFla = tfs->make<TH1F>("hNTrackperFla","hNTrackperFla",30,-5,25);
  hNTraperFla->GetXaxis()->SetTitle("N^{o} of CRTTrack per Flash");
  hNTraperFla->GetYaxis()->SetTitle("Entries/bin");

  hTra_tl_len = tfs->make<TH2F>("hTra_tl_len","hTra_tl_len",120, 0, 1200, 120, 0, 120);
  hTra_tl_len->GetXaxis()->SetTitle("Track lenght (cm)");
  hTra_tl_len->GetYaxis()->SetTitle("Track time (ns)");
  hTra_tl_len->SetOption("COLZ"); 

  
  hZdiff = tfs->make<TH1F>("hZdiff","hZdiff",100,-500,500);
  hZdiff->GetXaxis()->SetTitle("ZTrack - ZFlash (cm)");
  hZdiff->GetYaxis()->SetTitle("Entries/bin");

  hYdiff = tfs->make<TH1F>("hYdiff","hYdiff",100,-500,500);
  hYdiff->GetXaxis()->SetTitle("YTrack - YFlash (cm)");
  hYdiff->GetYaxis()->SetTitle("Entries/bin");
 
  double inch =2.54; //inch in cm
  hBot = tfs->make<TH2F>("hBottom","Bottom",125,-700+205*inch,-700+205*inch+125*10.89,60,-300+50.4*inch,-300+50.4*inch+60*10.89);
  hBot->GetXaxis()->SetTitle("Lenght along the beam (cm)");
  hBot->GetYaxis()->SetTitle("Lenght along the drift (cm)");
  hBot->GetZaxis()->SetTitle("Entries/bin");
  hBot->SetOption("COLZ");

  hFT = tfs->make<TH2F>("hFeedthroughSide","Feedthrough Side",125,-704+205*inch,-704+205*inch+125*10.89,60,-308-19.1*inch,-308-19.1*inch+60*10.89);
  hFT->GetXaxis()->SetTitle("Lenght along the beam (cm)");
  hFT->GetYaxis()->SetTitle("Height (cm)");
  hFT->GetZaxis()->SetTitle("Entries/bin");
  hFT->SetOption("COLZ");

  hPipe = tfs->make<TH2F>("hPipeSide","Pipe Side",125,-704+205*inch,-704+205*inch+125*10.89,60,-294-19.1*inch,-294-19.1*inch+60*10.89);
  hPipe->GetXaxis()->SetTitle("Lenght along the beam (cm)");
  hPipe->GetYaxis()->SetTitle("Height (cm)");
  hPipe->GetZaxis()->SetTitle("Entries/bin");
  hPipe->SetOption("COLZ");

  hTop = tfs->make<TH2F>("hTop","Top",125,-701+205*inch,-701+205*inch+125*11.38,80,2-170-300+50.4*inch,2-170-300+50.4*inch+80*11.38);
  hTop->GetXaxis()->SetTitle("Lenght along the beam (cm)");
  hTop->GetYaxis()->SetTitle("Lenght along the drift (cm)"); 
  hTop->GetZaxis()->SetTitle("Entries/bin"); 
  hTop->SetOption("COLZ");


  hTFvsTT = tfs->make<TH1F>("hTFvsTT","hTFvsTT",1000000,0,10000000);//1ms max
  hTFvsTT->GetXaxis()->SetTitle("Track time - Flash time (ns)");
  hTFvsTT->GetYaxis()->SetTitle("Entries/bin");

  hMulFT = tfs->make<TH1F>("hMulFT","hMulFT",50,0,50);//
  hMulFT->GetXaxis()->SetTitle("Multiplicity (Tracks per Flash)");
  hMulFT->GetYaxis()->SetTitle("Entries/bin");


  hMulFTvsTdis = tfs->make<TH2F>("hMulFTvsTdis","hMulFTvsTdis",50,0,50,1000000,0,10000000);
  hMulFTvsTdis->GetXaxis()->SetTitle("Multiplicity (Tracks per Flash)");
  hMulFTvsTdis->GetYaxis()->SetTitle("Track time - Flash time (ns)");
  hMulFTvsTdis->GetZaxis()->SetTitle("Entries/bin");
  hMulFTvsTdis->SetOption("COLZ");
  */
}

void T0recoCRTHit::endJob()
{
  // Implementation of optional member function here.
  
  
  /*	  
  //OLD
  //uint32_t Hit_sec = hitlist[ah]->ts0_s;
  //uint32_t Flash_sec = evt_timeGPS_sec;
  
  uint32_t Hit_nsec = hitlist[ah]->ts1_ns + fHardDelay_;
  //uint32_t Flash_nsec = Timeflash * 1000;
  
  int dif_sec = Flash_sec - Hit_sec;
  int dif_nsec = Flash_nsec - Hit_nsec;
  int dif_secABS = std::abs(dif_sec);
  int dif_nsecABS = std::abs(dif_nsec);
  //OLD
  
  if( (dif_secABS<3)  &&  (dif_nsecABS<1000 )  ){//E
  
  hTFvsTH->Fill(dif_nsec);
  hTFvsTH_2d->Fill(dif_sec , dif_nsecABS);
  
  hTFvsTH_t0->Fill(Timeflash_ns_GPS - hitlist[ah]->ts0_ns);
  hTFvsTH_t0_2d->Fill(Timeflash_ns_GPS - hitlist[ah]->ts0_ns, dif_nsec);
  hTFvsTH_t0_plane->Fill(hitlist[ah]->plane, Timeflash_ns_GPS - hitlist[ah]->ts0_ns);
  
  if(fverbose==1){
  std::cout<<"Flash_sec - Hit_sec: "<<Flash_sec - Hit_sec<<std::endl;
  std::cout<<"Flash_nsec - Hit_nsec: "<<Flash_nsec - Hit_nsec<<std::endl;
  std::cout<<"Flash_nsec - Hit_nsec: "<<dif_secABS<<std::endl;
  getchar();
  }
  
  }//E
  //OLD
  
  */
  
  
}

DEFINE_ART_MODULE(T0recoCRTHit)


