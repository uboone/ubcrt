////////////////////////////////////////////////////////////////////////////
/// Class:       CRTSimHitProducer
/// Module Type: producer
/// File:        CRTSimHitProducer_module.cc
///
/// Author:         Thomas Brooks
/// E-mail address: tbrooks@fnal.gov
///
/// Modified from CRTSimHitProducer by Thomas Warburton.
/////////////////////////////////////////////////////////////////////////////

#include "ubobj/CRT/CRTSimData.hh"
#include "ubobj/CRT/CRTHit.hh"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <map>
#include <iterator>
#include <algorithm>
#include <vector>

// LArSoft
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TGeoManager.h"

namespace {
  // Local namespace for local functions
    //variables to look up plane index, feb index and strip orientation for modules
    const short mod2plane[73] = {
      0,0,0,0,0,0,0,0,0,1,  //0-9
      1,1,1,1,1,1,1,1,1,1,  //10-19
      1,1,2,2,2,2,2,2,2,2,  //20-29
      2,2,2,2,2,2,2,2,2,2,  //30-39
      2,2,2,2,2,2,2,2,2,3,  //40-49
      3,3,3,3,3,3,3,3,3,3,  //50-59  
      3,3,3,3,3,3,3,3,3,3,  //60-69  
      3,3,3};               //70-72  
    const short mod2feb[73] = {
      24,23,22,17,14,18,19,12,11,52,  //0-9
      31,29,28,27,26,30,61,59,57,60,                  //10-19
      58,56,32,38,36,35,34,33,37,45,                  //20-29
      44,43,42,41,40,39,55,54,53,51,                  //30-39
      49,47,21,16,50,48,46,20,15,107,                 //40-49
      106,105,109,108,112,111,195,123,124,125,        //50-59
      126,129,115,114,113,116,119,121,127,128,        //60-69
      117,120,118};                                   //70-72
    const short mod2orient[73] = {  // 0=x, 1=y, 2=z
      2,2,2,0,0,0,0,2,0,1,  //0-9
      1,1,1,1,1,1,2,2,2,2,  //10-19
      2,2,1,1,1,1,1,1,1,1,  //20-29
      1,1,1,1,1,1,2,2,2,2,  //30-39
      2,2,2,2,2,2,2,2,2,0,  //40-49
      0,0,0,0,0,0,0,0,0,0,  //50-59
      0,2,2,2,2,2,2,2,2,0,  //60-69
      2,2,2};               //70-72

    /*  not currently used
    const short mod2channelorder[73] = {  //=0 if MCstrip ordering is increasing with coordinate, =1 if decreasing
      0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,
      0,0,0,1,1,1,1,1,1,1,
      1,0,1,1,1,1,1,1,1,0,
      1,1,1};
    const short mod2channelflip[73] = {  //=0 if strip ordering is the same for data and mc, =1 if opposite
      0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,1,
      1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,0,
      0,0,0,1,1,1,1,1,1,1,
      1,0,1,1,1,1,1,1,1,0,
      1,1,1};
    */


    
}
namespace crt{

  struct CRTStrip {
    double t0; 
    uint32_t channel; 
    double x; 
    double ex; 
    int id1; 
    int id2; 
    double pes1;
    double pes2;
    int module;
    int stripch;
  };
  
  class CRTSimHitProducer : public art::EDProducer {
  public:

    explicit CRTSimHitProducer(fhicl::ParameterSet const & p);

    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTSimHitProducer(CRTSimHitProducer const &) = delete;
    CRTSimHitProducer(CRTSimHitProducer &&) = delete;
    CRTSimHitProducer & operator = (CRTSimHitProducer const &) = delete; 
    CRTSimHitProducer & operator = (CRTSimHitProducer &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;

    void endJob() override;

    void reconfigure(fhicl::ParameterSet const & p);

    std::vector<double> ChannelToLimits(CRTStrip strip);

    std::vector<double> CrtOverlap(std::vector<double> strip1, std::vector<double> strip2);


    crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t,std::vector<std::pair<int,float>>> tpesmap, 
			   float peshit, double time, int plane,
			   double x, double ex, double y, double ey, double z, double ez); 

  private:

    // Params got from fcl file.......
    art::InputTag fCrtModuleLabel;      ///< name of crt producer
    bool          fVerbose;             ///< print info
    double        fTimeCoincidenceLimit;///< minimum time between two overlapping hit crt strips [ticks]
    double        fQPed;                ///< Pedestal offset of SiPMs [ADC]
    double        fQSlope;              ///< Pedestal slope of SiPMs [ADC/photon]
    bool          fUseReadoutWindow;    ///< Only reconstruct hits within readout window
    bool          fRequireStripOverlap;    ///< Make hits only if strips overlap
    double        fCRTClockFreq; // in GHz
    bool fTransAttenCorr;     // correct pes and position for transverse attenuation
    bool fLongAttenCorr;     // correct pes and time for longitudinal attenuation
   
    // Other variables shared between different methods.
    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    art::ServiceHandle<geo::AuxDetGeometry> fAuxDetGeoService;
    const geo::AuxDetGeometry* fAuxDetGeo;
    const geo::AuxDetGeometryCore* fAuxDetGeoCore;

  }; // class CRTSimHitProducer
    
  CRTSimHitProducer::CRTSimHitProducer(fhicl::ParameterSet const & p)
  // Initialize member data here, if know don't want to reconfigure on the fly
  {
    // Call appropriate produces<>() functions here.
    produces< std::vector<crt::CRTHit> >();
    
    // Get a pointer to the geometry service provider
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fAuxDetGeo = &(*fAuxDetGeoService);
    fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();

    reconfigure(p);

  } // CRTSimHitProducer()

  void CRTSimHitProducer::reconfigure(fhicl::ParameterSet const & p)
  {
    fCrtModuleLabel       = (p.get<art::InputTag> ("CrtModuleLabel")); 
    fVerbose              = (p.get<bool> ("Verbose"));
    fTimeCoincidenceLimit = (p.get<double> ("TimeCoincidenceLimit"));
    fQPed                 = (p.get<double> ("QPed"));
    fQSlope               = (p.get<double> ("QSlope"));
    fUseReadoutWindow     = (p.get<bool> ("UseReadoutWindow"),false);
    fRequireStripOverlap  = (p.get<bool> ("RequireStripOverlap"),false);
    fCRTClockFreq         = (p.get<double> ("CRTClockFreq"),1.0);  
    fTransAttenCorr       = (p.get<bool> ("TransAttenCorr"),true);
    fLongAttenCorr        = (p.get<bool> ("LongAttenCorr"),true);

  }

  void CRTSimHitProducer::beginJob()
    {
    if(fVerbose){std::cout<<"----------------- CRT Hit Reco Module -------------------"<<std::endl;}
    
  } // beginJob()
    
  void CRTSimHitProducer::produce(art::Event & event)
  {

    // const geo::AuxDetGeometry* geometry = &*fAuxDetGeoService;
    // const geo::AuxDetGeometryCore* geoServiceProvider = geometry->GetProviderPtr();

    int nHits = 0;

    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    // Detector properties
    double readoutWindow  = (double)fDetectorProperties->ReadOutWindowSize();
    double driftTimeTicks = 2.0*(2.*fGeometryService->DetHalfWidth()+3.)/fDetectorProperties->DriftVelocity();

    // Retrieve list of CRT hits
    art::Handle< std::vector<crt::CRTSimData>> crtListHandle;
    std::vector<art::Ptr<crt::CRTSimData> > crtList;
    if (event.getByLabel(fCrtModuleLabel, crtListHandle))
      art::fill_ptr_vector(crtList, crtListHandle);

    // Place to store CRThits as they are created
    std::unique_ptr<std::vector<crt::CRTHit>> CRTHitcol( new std::vector<crt::CRTHit>);
    std::map<int, std::vector<CRTStrip>> taggerStrips;

    // Loop over all the SiPM hits in 2 (should be in pairs due to trigger)
    if(fVerbose) std::cout<<"Number of SiPM hits = "<<crtList.size()<< std::endl;

    for (size_t i = 0; i < crtList.size(); i+=2){


      art::Ptr<crt::CRTSimData> thisSiPM1 = crtList[i];      
      art::Ptr<crt::CRTSimData> thisSiPM2 = crtList[i+1];      

      // Get the time, channel, center and width
      double t1 = (double)((int)(thisSiPM1->fT0))*fCRTClockFreq*1e-3;
      if(fUseReadoutWindow){
        if(!(t1 >= -driftTimeTicks && t1 <= readoutWindow)) continue;
      }
      uint32_t channel = thisSiPM1->fChannel/2;
      // channel here is really the AuxDetID
      std::string name = fGeometryService->AuxDet(channel).TotalVolume()->GetName();
      int strip=0; int module=0;
      sscanf(name.c_str(),"volAuxDet_Module_%d_strip_%d",&module,&strip);
      if (fVerbose) std::cout << " channel " << channel << " strip " << strip << 
		      " module " << module << " name " << name << std::endl;
      if (module<0 || module>72) {
	std::cout << " module out of range " << module << std::endl;
	module=0;
      }

      // width of all strips is 10.8 cm except the strip in the top plane modules, 
      //         which are 11.2 cm wide
      //  . . . . hardcoding is ugly.
      int thisplane = mod2plane[module];
      double width = 10.8;  
      if (thisplane==3) width = 11.2;
      int this_orient=0;
      if (mod2orient[module]==2) this_orient=1;
      int thistagger=2*thisplane+this_orient;
      if (fVerbose) std::cout << "plane " << thisplane << " orient " << this_orient << " tagger " << 
		      thistagger << std::endl;

      // track IDs from geant
      int id1 = thisSiPM1->fTrackID;
      int id2 = thisSiPM2->fTrackID; 

      // Get the time of hit on the second SiPM in us
      double t2 = (double)(int)thisSiPM2->fT0*fCRTClockFreq*1e-3;
      // Calculate the number of photoelectrons at each SiPM
      double npe1 = ((double)thisSiPM1->fADC - fQPed)/fQSlope;
      double npe2 = ((double)thisSiPM2->fADC - fQPed)/fQSlope;
      // Calculate the hit position across the width of the strip (between two sipms)
      double x,ex;
      if (fTransAttenCorr) {
	x = (width/2.)*atan(log(1.*npe2/npe1)) + (width/2.);
	double normx = x + 0.344677*x - 1.92045;
	ex = 1.92380e+00+1.47186e-02*normx-5.29446e-03*normx*normx;
      }
      else {x = -0.5*width+ (width)*(npe2/(npe1+npe2)); ex=1.0;}

      if(fVerbose) std::cout << "strip hit: time (us) "<< t1 << " channel " <<  channel <<
      		     " tran pos " << x <<  "  id1 " << id1 << " id2 " <<
      		     id2 << " total pe " <<  npe1+npe2  << std::endl;

      double ttime = (t1 + t2)/2.;

      CRTStrip stripHit = {ttime, channel, x, ex, id1, id2, npe1, npe2, module, strip};

      if (taggerStrips.count(thistagger)>0) {
	auto search = taggerStrips.find(thistagger);
	(search->second).push_back(stripHit);
      }
      else {
	std::vector<CRTStrip> newshvec; newshvec.push_back(stripHit);	
	taggerStrips.emplace(thistagger,newshvec);
      }


    }// end loop over sipm signals
    
    
    // Remove any duplicate (same channel and time) hit strips
    for(auto &tagStrip : taggerStrips){
      std::sort(tagStrip.second.begin(), tagStrip.second.end(),
                [](const CRTStrip & a, const CRTStrip & b) -> bool{
                  return (a.t0 < b.t0) || 
		    ((a.t0 == b.t0) && (a.channel < b.channel));
                });
      // Remove hits with the same time and channel
      tagStrip.second.erase(std::unique(tagStrip.second.begin(), tagStrip.second.end(),
					[](const CRTStrip & a, const CRTStrip & b) -> bool{
					  return a.t0 == b.t0 && a.channel == b.channel;
					}), tagStrip.second.end());
    }// end loop over tagger lists of strip hits.
    
    if (fVerbose) std::cout << " lists of strip hits complete "  <<
		    taggerStrips.count(0) << " " <<taggerStrips.count(1) << " " <<
		    taggerStrips.count(2) << " " <<taggerStrips.count(3) << " " <<
		    taggerStrips.count(4) << " " <<taggerStrips.count(5) << " " <<
		    taggerStrips.count(6) << " " <<taggerStrips.count(7) << std::endl;

    //loop over planes, look for hits at the same time from strips of opposite orientation.  If the
    //   strips overlap, make a CRT hit.
    for (int ip=0;ip<4;++ip) {   
      int index1=2*ip; int index2=2*ip+1;
      if (taggerStrips.count(index1)>0 && taggerStrips.count(index2)>0) {
	auto search = taggerStrips.find(index1);
	std::vector<CRTStrip> tpo1=search->second;
	auto search2 = taggerStrips.find(index2);
	std::vector<CRTStrip> tpo2=search2->second;
	for (size_t hit_i = 0; hit_i < tpo1.size(); hit_i++){
	  CRTStrip thisstrip1 = tpo1[hit_i];	
	  int thismodule = thisstrip1.module;
	  int thisplane = mod2plane[thismodule];
	  int stripdir = mod2orient[thismodule];
	  double t0_1 = thisstrip1.t0;
	  std::vector<double> limits1 = ChannelToLimits(thisstrip1);
	  // Loop over hits in orthogonal strips
	  for (size_t hit_j = 0; hit_j < tpo2.size(); hit_j++){
	    //check for opposite strip orientation
	    CRTStrip thisstrip2 = tpo2[hit_j];	
	    int thismodule2 = thisstrip2.module;
	    //	    int thisplane2 = mod2plane[thismodule2];
	    //	    int stripdir2 = mod2orient[thismodule2];
	    double t0_2 = thisstrip2.t0;
	    if (std::abs(t0_1 - t0_2)<fTimeCoincidenceLimit) {
	      // Average the time
	      double time = (t0_1 + t0_2)/2;
	      //get FEB IDs
	      std::vector<uint8_t> tfeb_id; 
	      tfeb_id.push_back(mod2feb[thismodule]);
	      tfeb_id.push_back(mod2feb[thismodule2]);
	      std::vector<std::pair<int,float>> myvec1,myvec2;
	      myvec1.push_back(std::pair<int,float>(thisstrip1.stripch*2,thisstrip1.pes1));
	      myvec1.push_back(std::pair<int,float>(thisstrip1.stripch*2+1,thisstrip1.pes2));
	      myvec2.push_back(std::pair<int,float>(thisstrip2.stripch*2,thisstrip2.pes1));
	      myvec2.push_back(std::pair<int,float>(thisstrip2.stripch*2+1,thisstrip2.pes2));
	      std::map<uint8_t, std::vector<std::pair<int,float>>> mymap;	      
	      uint8_t if1 = mod2feb[thismodule];
	      mymap.insert(std::pair<uint8_t, std::vector<std::pair<int,float>>>(if1,myvec1));
	      uint8_t if2 = mod2feb[thismodule2];
	      mymap.insert(std::pair<uint8_t, std::vector<std::pair<int,float>>>(if2,myvec2));
	      double petot = thisstrip1.pes1+ thisstrip1.pes2 + thisstrip2.pes1+ thisstrip2.pes2;

	      //check for strip overlap
	      if (fVerbose) std::cout << " checking overlap " << hit_i << " " << hit_j  << std::endl;
	      std::vector<double> limits2 = ChannelToLimits(thisstrip2);
	      std::vector<double> overlap = CrtOverlap(limits1, limits2);
	      if (overlap[0] != -99999)  {
		// Calculate the mean and error in x, y, z
		TVector3 mean((overlap[0] + overlap[1])/2., 
			      (overlap[2] + overlap[3])/2., 
			      (overlap[4] + overlap[5])/2.);
		TVector3 error(std::abs((overlap[1] - overlap[0])/2.), 
			       std::abs((overlap[3] - overlap[2])/2.), 
			       std::abs((overlap[5] - overlap[4])/2.));
		if (thisplane==0 || thisplane==3) error.SetY(1.0);
		else error.SetX(1.0);
		// Create a CRT hit
		crt::CRTHit crtHit = FillCrtHit(tfeb_id, mymap, petot, time, thisplane, mean.X(),
						error.X(), mean.Y(), error.Y(), mean.Z(), error.Z());
		CRTHitcol->push_back(crtHit);
		nHits++;
		if (fVerbose) std::cout << "hit created: time " << time << " x " <<  mean.X() << 
				" y " << mean.Y() << " z " <<  mean.Z() << std::endl;
		
	      }
	      else if (!fRequireStripOverlap) {	     
		TVector3 mean,error;
		if (thisplane==0 || thisplane==3) { // top or bot planes at constant y
		  if (stripdir==0) {  //strip1 length along x
		    mean.SetZ(0.5*(limits1[4]+limits1[5]));
		    error.SetZ(0.5*std::abs(limits1[5]-limits1[4]));
		    mean.SetX(0.5*(limits2[0]+limits2[1]));
		    error.SetX(0.5*std::abs(limits2[1]-limits2[0]));
		    mean.SetY(0.25*(limits1[2]+limits1[3]+limits2[2]+limits2[3]));
		    error.SetY(0.5);
		  }
		  else {
		    mean.SetZ(0.5*(limits2[4]+limits2[5]));
		    error.SetZ(0.5*std::abs(limits2[5]-limits2[4]));
		    mean.SetX(0.5*(limits1[0]+limits1[1]));
		    error.SetX(0.5*std::abs(limits1[1]-limits1[0]));
		    mean.SetY(0.25*(limits1[2]+limits1[3]+limits2[2]+limits2[3]));
		    error.SetY(0.5);
		  }
		}
		else {  // side planes at constant x
		  if (stripdir==1) { //strip1 length along y
		    mean.SetZ(0.5*(limits1[4]+limits1[5]));
		    error.SetZ(0.5*std::abs(limits1[5]-limits1[4]));
		    mean.SetY(0.5*(limits2[3]+limits2[2]));
		    error.SetY(0.5*std::abs(limits2[3]-limits2[2]));
		    mean.SetX(0.25*(limits1[0]+limits1[1]+limits2[0]+limits2[1]));
		    error.SetX(0.5);
		  }
		  else {
		    mean.SetZ(0.5*(limits2[4]+limits2[5]));
		    error.SetZ(0.5*std::abs(limits2[5]-limits2[4]));
		    mean.SetY(0.5*(limits1[3]+limits1[2]));
		    error.SetY(0.5*std::abs(limits1[3]-limits1[2]));
		    mean.SetX(0.25*(limits1[0]+limits1[1]+limits2[0]+limits2[1]));
		    error.SetX(0.5);
		  }
		}
		// Create a CRT hit
		crt::CRTHit crtHit = FillCrtHit(tfeb_id, mymap, petot, time, thisplane, mean.X(), 
						error.X(), mean.Y(), error.Y(), mean.Z(), error.Z());
		CRTHitcol->push_back(crtHit);
		nHits++;
		if (fVerbose) std::cout << "hit created: time " << time << " x " <<  mean.X() << 
				" y " << mean.Y() << " z " <<  mean.Z() << std::endl;
	      }// end create hit if strip overlap requirement is met (or not set)
	      }// if coincidence in time
	      }  //end loop over pairing
	    }  // end loop over strip hits in that sub-plane
	  }// end if there are hits
	}  // end loop over planes
	
	event.put(std::move(CRTHitcol));
	
	if(fVerbose) std::cout<<"Number of CRT hits produced = "<<nHits<<std::endl;
	
      } // produce()
    
      void CRTSimHitProducer::endJob()
      {
    
      }

      // Function to calculate the strip position limits in real space from channel
      std::vector<double> CRTSimHitProducer::ChannelToLimits(CRTStrip strHit)  {

    
	// uncommenting this causes it to seg fault, not sure why.
	// const geo::AuxDetGeo stripGeoL = fGeometryService->AuxDet(strHit.channel);
	// double mycenter[3];     stripGeoL.GetCenter(mycenter);
	// if (fVerbose) std::cout << "center " << mycenter[0] << " " << mycenter[1] << " " << 
	// 		    mycenter[2]<< std::endl;

	const geo::AuxDetSensitiveGeo stripGeoL = (fGeometryService->AuxDet(strHit.channel)).SensitiveVolume(0);
	// double mycenter[3];     stripGeoL.GetCenter(mycenter);
	// if (fVerbose) std::cout << "center " << mycenter[0] << " " << mycenter[1] << " " << 
	// 		    mycenter[2]<< std::endl;


	double thisx1,thisy1,thisz1,thisx2,thisy2,thisz2;

	double halfheight = stripGeoL.HalfHeight();
	double halfwidth = stripGeoL.HalfWidth1();
	double halflength = 0.5*(stripGeoL.Length());

	int stripaxis = mod2orient[strHit.module];    // 0=x, 1=y, 2=z
	int stripplane = mod2plane[strHit.module];    // 0=bot, 1=anode, 2=cathode, 3=top

	if (fVerbose) {
	  std::cout << " plane " << stripplane << " axis " << stripaxis 
		    << " width " << halfwidth<< " height " << 
	    halfheight  << " length " << halflength << std::endl;
	  std::cout << strHit.x << " " << strHit.ex << std::endl;
	}
	if (stripaxis==0) {
	  thisx1= halfwidth;
	  thisy1= halfheight;
	  thisz1= -1.0*halflength+strHit.x-strHit.ex;
	  thisx2= -1.0*halfwidth;
	  thisy2= -1.0*halfheight;
	  thisz2= -1.0*halflength+strHit.x+strHit.ex;
	}
	else if (stripaxis==1) {
	  thisx1= halfwidth;
	  thisy1= halfheight;
	  thisz1= -1.0*halflength+strHit.x-strHit.ex;
	  thisx2= -1.0*halfwidth;
	  thisy2= -1.0*halfheight;
	  thisz2= -1.0*halflength+strHit.x+strHit.ex;
	}
	else {  // stripaxis is z
	  if (stripplane==0 || stripplane==3) {
	    thisx1= -1.0*halfwidth+strHit.x-strHit.ex;
	    thisy1= halfheight;
	    thisz1= halflength;
	    thisx2= -1.0*halfwidth+strHit.x+strHit.ex;
	    thisy2= -1.0*halfheight;
	    thisz2= -1.0*halflength;
	  }
	  else {
	    thisx1= halfwidth;
	    thisy1= -1.0*halfheight+strHit.x-strHit.ex;
	    thisz1= halflength;
	    thisx2= -1.0*halfwidth;
	    thisy2= -1.0*halfheight+strHit.x+strHit.ex;
	    thisz2= -1.0*halflength;
	  }
	}
	// give more room for two crossing strips to overlap in the axis perpendicular to the plane
	//  allowed difference will be +/- 10.0 cm now instead of +/- 0.5 cm
	// hardcoding the geometry is ugly
	if (stripplane==3 || stripplane==0) { thisy1+=10.0; thisy2-=10.0; }
	else {thisx1+=10.0; thisx2-=10.0;}

	double w1[3] = {0,0,0};
	double w2[3] = {0,0,0};
	const double l1[3] = {thisx1,thisy1,thisz1};
	stripGeoL.LocalToWorld(l1, w1);
	const double l2[3] = {thisx2,thisy2,thisz2};
	stripGeoL.LocalToWorld(l2, w2);
	// if (fVerbose) {
	//   std::cout << " l1 " << l1[0] << " " << l1[1] << " " << l1[2] << std::endl;
	//   std::cout << " w1 " << w1[0] << " " << w1[1] << " " << w1[2] << std::endl;
	//   std::cout << " l2 " << l2[0] << " " << l2[1] << " " << l2[2] << std::endl;
	//   std::cout << " w2 " << w2[0] << " " << w2[1] << " " << w2[2] << std::endl;
	// }

    // Use this to get the limits in the two variable directions
    std::vector<double> limits = {std::min(w1[0],w2[0]), std::max(w1[0],w2[0]), 
                                  std::min(w1[1],w2[1]), std::max(w1[1],w2[1]), 
                                  std::min(w1[2],w2[2]), std::max(w1[2],w2[2])};
    // if (fVerbose) std::cout << strHit.module << " limits " << limits[0] << " " << limits[1] << " " << limits[2] << " "
    // 			     << limits[3] << " " << limits[4] << " " << limits[5] << " " << std::endl;

    return(limits);
  } // ChannelToLimits

		      
  // Function to calculate the overlap between two crt strips
    std::vector<double> CRTSimHitProducer::CrtOverlap(std::vector<double> strip1,std::vector<double> strip2) {
     
    double minX = std::max(strip1[0], strip2[0]);
    double maxX = std::min(strip1[1], strip2[1]);
    double minY = std::max(strip1[2], strip2[2]);
    double maxY = std::min(strip1[3], strip2[3]);
    double minZ = std::max(strip1[4], strip2[4]);
    double maxZ = std::min(strip1[5], strip2[5]);

    std::vector<double> null = {-99999, -99999, -99999, -99999, -99999, -99999};
    std::vector<double> overlap = {minX, maxX, minY, maxY, minZ, maxZ};
    if ((minX<maxX && minY<maxY) || (minX<maxX && minZ<maxZ) || (minY<maxY && minZ<maxZ)) return(overlap);

    return(null);

  } 



    crt::CRTHit CRTSimHitProducer::FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap, float peshit,double time, int plane, double x, double ex, double y, double ey, double z, double ez){

    crt::CRTHit crtHit;
    crtHit.feb_id = tfeb_id;
    crtHit.pesmap = tpesmap;
    crtHit.peshit = peshit;
    crtHit.ts0_s_corr = 0;
    crtHit.ts0_ns = 0;
    crtHit.ts0_ns_corr = 0;
    crtHit.ts1_ns = time * 0.5 * 10e3;
    crtHit.ts0_s = time * 0.5 * 10e-6; 
    crtHit.plane = plane;
    crtHit.x_pos = x;
    crtHit.x_err = ex;
    crtHit.y_pos = y; 
    crtHit.y_err = ey;
    crtHit.z_pos = z;
    crtHit.z_err = ez;
    return crtHit;
  }


  DEFINE_ART_MODULE(CRTSimHitProducer)

    }// namespace crt

namespace {


}
