////////////////////////////////////////////////////////////////////////////
/// Class:       CRTSimHitCorr
/// Module Type: producer
/// File:        CRTSimHitCorr_module.cc
///
/// Author:         Michelle Stancari
/// E-mail address: mstancar@fnal.gov
///
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
#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

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


  const short feb2mod[200] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //0-9
    -1,8,7,-1,4,48,43,3,5,6, //10-19
    47,42,2,1,0,-1,14,13,12,11,  //20-29
    15,10,22,27,26,25,24,28,23,35,  //30-39
    34,33,32,31,30,29,46,41,45,40,  //40-49
    44,39,9,38,37,36,21,18,20,17,   //50-59
    19,16,-1,-1,-1,-1,-1,-1,-1,-1, //60-69
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //70-79
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //80-89
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //90-99
    -1,-1,-1,-1,-1,51,50,49,53,52, //100-109
    -1,55,54,64,63,62,65,70,72,66, //110-119
    71,67,-1,57,58,59,60,68,69,61, //120-129
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //130-139
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //140-149
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //150-159
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //160-169
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //170-179
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //180-189
    -1,-1,-1,-1,-1,56,-1,-1,-1,-1}; //190-199

  const float sipm_pos[73] = {  // index is module number from gdml
    628.952, 628.952, 628.952, -134.384, -134.384, 390.616, 390.616, 515.952, -134.384, -218.514,  //0-9
    -218.514, -218.514, -218.514, -218.514, -218.514, -218.514, 312.7, 716.7, 1120.7, 312.7,                  //10-19
    716.7, 1120.7, 291.486, 291.486, 291.486, 291.486, 291.486, 291.486, 291.486, -218.514,                  //20-29
    -218.514, -218.514, -218.514, -218.514, -218.514, -218.514, 327.7, 724.2, 1120.7, 188.7,                  //30-39
    416.7, 644.7, 872.7, 1100.7, 188.7, 416.7, 644.7, 827.7, 1100.7, -230.0,                 //40-49
    -230.0, -230.0, -230.0, 490.0, 490.0, 490.0, 310.0, 310.0, 310.0, -50.0,        //50-59
    -50.0, -80.0, 1180.0, 1180.0, 1180.0, 820.0, 100.0, 100.0, 640.0, 130.0,        //60-69
    820.0, 100.0, 820.0};                                   //70-72

/*
    const short mod2feb[73] = {
      24,23,22,17,14,18,19,12,11,52,  //0-9
      31,29,28,27,26,30,61,59,57,60,                  //10-19
      58,56,32,38,36,35,34,33,37,45,                  //20-29
      44,43,42,41,40,39,55,54,53,51,                  //30-39
      49,47,21,16,50,48,46,20,15,107,                 //40-49
      106,105,109,108,112,111,195,123,124,125,        //50-59
      126,129,115,114,113,116,119,121,127,128,        //60-69
      117,120,118};                                   //70-72
*/
    const short mod2orient[73] = {  // 0=x, 1=y, 2=z
      2,2,2,0,0,0,0,2,0,1,  //0-9
      1,1,1,1,1,1,2,2,2,2,  //10-19
      2,2,1,1,1,1,1,1,1,1,  //20-29
      1,1,1,1,1,1,2,2,2,2,  //30-39
      2,2,2,2,2,2,2,2,2,0,  //40-49
      0,0,0,0,0,0,0,0,0,0,  //50-59
      0,2,2,2,2,2,2,2,2,0,  //60-69
      2,2,2};               //70-72

    const float mod2length[73] = { 
      346.0,346.0,346.0,259.6,259.6,259.6,259.6,227.0,227.0,346.0,  //0-9
      346.0,346.0,346.0,346.0,346.0,346.0,403.8,403.8,403.8,403.8,  //10-19
      403.8,403.8,259.6,259.6,259.6,259.6,259.6,259.6,259.6,259.6,  //20-29
      259.6,259.6,259.6,259.6,259.6,259.6,396.2,396.2,396.2,227.0,  //30-39
      227.0,227.0,227.0,227.0,227.0,227.0,227.0,227.0,227.0,360.0,  //40-49
      360.0,360.0,360.0,360.0,360.0,360.0,180.0,180.0,180.0,180.0,  //50-59
      180.0,180.0,365.0,365.0,365.0,365.0,365.0,365.0,180.0,180.0,  //60-69
      365.0,365.0,365.0};               //70-72

    const short mod2end[73] = {  // -1 means sipm is a higher coordinate value than the center of the strip
      -1,-1,-1,+1,+1,-1,-1,+1,+1,+1,  //0-9
      +1,+1,+1,+1,+1,+1,-1,-1,-1,-1,  //10-19
      -1,-1,-1,-1,-1,-1,-1,-1,-1,+1,  //20-29
      +1,+1,+1,+1,+1,+1,-1,-1,-1,-1,  //30-39
      -1,-1,-1,-1,-1,-1,-1,-1,-1,+1,  //40-49
      +1,+1,+1,-1,-1,-1,-1,-1,-1,+1,  //50-59
      +1,+1,-1,-1,-1,-1,+1,+1,-1,+1,  //60-69
      -1,+1,-1};               //70-72

namespace crt{

  
  class CRTSimHitCorr : public art::EDProducer {
  public:

    explicit CRTSimHitCorr(fhicl::ParameterSet const & p);

    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTSimHitCorr(CRTSimHitCorr const &) = delete;
    CRTSimHitCorr(CRTSimHitCorr &&) = delete;
    CRTSimHitCorr & operator = (CRTSimHitCorr const &) = delete; 
    CRTSimHitCorr & operator = (CRTSimHitCorr &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;

    void endJob() override;

    void reconfigure(fhicl::ParameterSet const & p);


    crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t,std::vector<std::pair<int,float>>> tpesmap, 
			   float peshit, double time1, double time2, double time3, double time4, double time5, int plane,
			   double x, double ex, double y, double ey, double z, double ez); 

  private:

    // Params got from fcl file.......
    art::InputTag fCrtHitsIn_Label;      ///< name of crt producer
    bool          fScaleMCtime;             ///< turns off/on bug fix for hit times
    float fHitThreshold;
    float fStripThreshold;
    float fSiPMThreshold;
    float fPEscaleFactor;
    bool fRemoveBottomHits;
    bool fApplyDetectorResponse;
    bool          fVerbose;             ///< print info
    CLHEP::HepRandomEngine& fEngine;
   

  }; // class CRTSimHitCorr
    
  CRTSimHitCorr::CRTSimHitCorr(fhicl::ParameterSet const & p)
    : EDProducer{p},     
    fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "crt", p, "Seed"))

  // Initialize member data here, if know don't want to reconfigure on the fly
  {
    // Call appropriate produces<>() functions here.
    produces< std::vector<crt::CRTHit> >();
    
    // fEngine = art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "crt", p, "Seed");

    reconfigure(p);

  } // CRTSimHitCorr()

  void CRTSimHitCorr::reconfigure(fhicl::ParameterSet const & p)
  {
    fCrtHitsIn_Label       = (p.get<art::InputTag> ("CrtHitsIn_Label","crtsimhit")); 
    fScaleMCtime           = (p.get<bool>("ScaleMCtime",false));
    fHitThreshold           = (p.get<float>("HitThreshold",0.0));
    fStripThreshold           = (p.get<float>("StripThreshold",0.0));
    fSiPMThreshold           = (p.get<float>("SiPMThreshold",0.0));
    fPEscaleFactor           = (p.get<float>("PEscaleFactor",-1.0));
    fRemoveBottomHits           = (p.get<bool>("RemoveBottomHits",true));
    fApplyDetectorResponse           = (p.get<bool>("ApplyDetectorResponse",false));
    fVerbose              = (p.get<bool> ("Verbose",false));
    //

  }

  void CRTSimHitCorr::beginJob()
    {
    if(fVerbose){std::cout<<"----------------- CRT Hit Reco Module -------------------"<<std::endl;}
    
  } // beginJob()
    
  void CRTSimHitCorr::produce(art::Event & event)
  {

    int nHits = 0;

    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    // Place to store corrected CRThits as they are created
    std::unique_ptr<std::vector<crt::CRTHit>> CRTHitOutCol( new std::vector<crt::CRTHit>);
    // Retrieve list of CRT hits
    art::Handle< std::vector<crt::CRTHit>> crtHitsInHandle;
    event.getByLabel(fCrtHitsIn_Label, crtHitsInHandle);
   //check to make sure the data we asked for is valid
    if(!crtHitsInHandle.isValid()){
      std::cout << "Run " << event.run() << ", subrun " << event.subRun()
		<< ", event " << event.event() << " has zero"
		<< " CRTHits " << " in module " << fCrtHitsIn_Label << std::endl;
      std::cout << std::endl;
      //add protection here
      event.put(std::move(CRTHitOutCol));
      return;
    }
    std::vector<crt::CRTHit> const& crtHitInList(*crtHitsInHandle);
    if(fVerbose) std::cout<<"Number of CRT hits read in= "<<crtHitInList.size()<< std::endl;



    for (size_t i = 0; i < crtHitInList.size(); i++){

      crt::CRTHit thisCrtHit = crtHitInList[i];

      std::vector<uint8_t> tfeb_id = thisCrtHit.feb_id; 
      double time1 = thisCrtHit.ts0_s;
      double time2 = thisCrtHit.ts0_s_corr;
      double time3 = thisCrtHit.ts0_ns;
      double time4 = thisCrtHit.ts0_ns_corr;
      double time5 = thisCrtHit.ts1_ns;
      
      int plane = thisCrtHit.plane;
      double x = thisCrtHit.x_pos;
      double ex = thisCrtHit.x_err;
      double y = thisCrtHit.y_pos;
      double ey = thisCrtHit.y_err;
      double z = thisCrtHit.z_pos;
      double ez= thisCrtHit.z_err;
      
      std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap=thisCrtHit.pesmap;
      float pestot = thisCrtHit.peshit;      

      int iKeepMe = 1;
      // only change/remove MC hits
      std::vector<std::pair<int,float>> test = tpesmap.find(tfeb_id[0])->second; 
      if (test.size()==2)  { // this is simulation
	//  apply corrections to fix bug in Oct2018 simulation
	if (fScaleMCtime) {
	  time5/=5;
	  time1/=5;
	}
	if (fApplyDetectorResponse) {	  
	  if (fPEscaleFactor<0) fPEscaleFactor=1.0; 
	  tpesmap.clear();
	  // first strip
	  double distToReadout=0;
	  int this_feb = tfeb_id[0];
	  int this_mod = feb2mod[this_feb];
	  if (this_mod<0 || this_mod>72) 
	    std::cout << "bad module number for feb " << this_feb << std::endl;
	  else {	    
	    if (mod2orient[this_mod]==0) distToReadout=mod2end[this_mod]*(x-sipm_pos[this_mod]);
	    else if (mod2orient[this_mod]==1) distToReadout=mod2end[this_mod]*(y-sipm_pos[this_mod]);
	    else distToReadout=mod2end[this_mod]*(z-sipm_pos[this_mod]);
	    if (distToReadout>mod2length[this_mod]) {
	      //	      std::cout << "  HMMM  module/feb " << this_mod << " " << this_feb << " " << distToReadout << " " << 
	      //	mod2length[this_mod] << std::endl;
	      distToReadout=mod2length[this_mod];
	    }
	    else if (distToReadout<0) distToReadout=0.0;
	  }	  
	  double b = 1085.0;
	  double pe_sf_A = b*b/pow(distToReadout+b,2.0);

	  std::map<uint8_t, std::vector<std::pair<int,float>>> tempmap=thisCrtHit.pesmap;
	  std::vector<std::pair<int,float>> pesA = tempmap.find(tfeb_id[0])->second; 
	  float pesA_0=(fPEscaleFactor*pe_sf_A)*pesA[0].second;
	  float pesB_0=(fPEscaleFactor*pe_sf_A)*pesA[1].second;
	  float pestot_0=pesA_0+pesB_0;
	  float distA=(pesA_0/pestot_0)*11.0;  // 11.0 is strip width in cm
	  float distB=(pesB_0/pestot_0)*11.0;  // 11.0 is strip width in cm
	  float absA=exp(-1.0*distA/8.5);
	  float absB=exp(-1.0*distB/8.5);
	  pesA_0=pestot_0*(absA)/(absA+absB);
	  pesB_0=pestot_0*(absB)/(absA+absB);
	  // smear these values
          float npe0 = CLHEP::RandPoisson::shoot(&fEngine, pesA_0);
          float npe1 = CLHEP::RandPoisson::shoot(&fEngine, pesB_0);
	  // SiPM and ADC response: Npe to ADC counts
          float pesA_sm = CLHEP::RandGauss::shoot(&fEngine, npe0, 0.085 * sqrt(npe0));
          float pesB_sm = CLHEP::RandGauss::shoot(&fEngine, npe1, 0.085 * sqrt(npe1));
	  std::pair<int,float> pesA0(pesA[0].first,pesA_sm);
	  std::pair<int,float> pesA1(pesA[1].first,pesB_sm);
	  std::vector<std::pair<int,float>> pesAnew;
	  pesAnew.push_back(pesA0); pesAnew.push_back(pesA1);
	  tpesmap.emplace(tfeb_id[0],pesAnew);
	  pestot = pesA_sm+pesB_sm;

	  //second strip
	  distToReadout=0;
	  this_feb = tfeb_id[1];
	  this_mod = feb2mod[this_feb];
	  if (this_mod<0 || this_mod>72) 
	    std::cout << "bad module number for feb " << this_feb << std::endl;
	  else {	    
	    if (mod2orient[this_mod]==0) distToReadout=mod2end[this_mod]*(x-sipm_pos[this_mod]);
	    else if (mod2orient[this_mod]==1) distToReadout=mod2end[this_mod]*(y-sipm_pos[this_mod]);
	    else distToReadout=mod2end[this_mod]*(z-sipm_pos[this_mod]);
	    if (distToReadout>mod2length[this_mod]) {
	      // std::cout << "  HMMM  module/feb " << this_mod << " " << this_feb << " " << distToReadout << " " << 
	      // 	mod2length[this_mod] << std::endl;
	      distToReadout=mod2length[this_mod];
	    }
	    else if (distToReadout<0) distToReadout=0.0;
	  }	  
	  double pe_sf_B = b*b/pow(distToReadout+b,2.0);

	  std::vector<std::pair<int,float>> pesB = tempmap.find(tfeb_id[1])->second; 
	  pesA_0=pe_sf_B*fPEscaleFactor*pesB[0].second;
	  pesB_0=pe_sf_B*fPEscaleFactor*pesB[1].second;
	  pestot_0=pesA_0+pesB_0;
	  distA=(pesA_0/pestot_0)*11.0;  // 11.0 is strip width in cm
	  distB=(pesB_0/pestot_0)*11.0;  // 11.0 is strip width in cm
	  absA=exp(-1.0*distA/8.5);
	  absB=exp(-1.0*distB/8.5);
	  pesA_0=pestot_0*(absA)/(absA+absB);
	  pesB_0=pestot_0*(absB)/(absA+absB);
	  // smear these values
          npe0 = CLHEP::RandPoisson::shoot(&fEngine, pesA_0);
          npe1 = CLHEP::RandPoisson::shoot(&fEngine, pesB_0);
	  // SiPM and ADC response: Npe to ADC counts
          pesA_sm = CLHEP::RandGauss::shoot(&fEngine, npe0, 0.085 * sqrt(npe0));
          pesB_sm = CLHEP::RandGauss::shoot(&fEngine, npe1, 0.085 * sqrt(npe1));
	  std::pair<int,float> pesB0(pesB[0].first,pesA_sm);
	  std::pair<int,float> pesB1(pesB[1].first,pesB_sm);
	  std::vector<std::pair<int,float>> pesBnew;
	  pesBnew.push_back(pesB0); pesBnew.push_back(pesB1);
	  tpesmap.emplace(tfeb_id[1],pesBnew);
	  pestot += pesA_sm+pesB_sm;

	} // if Apply Detector Response Flag is set
	else if (fPEscaleFactor>0) {
	  // scale total hit charge 
	  pestot*=fPEscaleFactor;
	  // more scaling of charge, maps are painful.
	  tpesmap.clear();
	  std::map<uint8_t, std::vector<std::pair<int,float>>> tempmap=thisCrtHit.pesmap;
	  std::vector<std::pair<int,float>> pesA = tempmap.find(tfeb_id[0])->second; 
	  std::pair<int,float> pesA0(pesA[0].first,fPEscaleFactor*pesA[0].second);
	  std::pair<int,float> pesA1(pesA[1].first,fPEscaleFactor*pesA[1].second);
	  std::vector<std::pair<int,float>> pesAnew;
	  pesAnew.push_back(pesA0); pesAnew.push_back(pesA1);
	  tpesmap.emplace(tfeb_id[0],pesAnew);
	  std::vector<std::pair<int,float>> pesB = tempmap.find(tfeb_id[1])->second; 
	  std::pair<int,float> pesB0(pesB[0].first,fPEscaleFactor*pesB[0].second);
	  std::pair<int,float> pesB1(pesB[1].first,fPEscaleFactor*pesB[1].second);
	  std::vector<std::pair<int,float>> pesBnew;
	  pesBnew.push_back(pesB0); pesBnew.push_back(pesB1);
	  tpesmap.emplace(tfeb_id[1],pesBnew);
	} // if scale factor >0 and det response turned off

	
	// remove bottom hits from FEB combinations not allowed in data
	if (fRemoveBottomHits) {
	  if ((tfeb_id[0]==11 && tfeb_id[1]!=12) || (tfeb_id[0]==12 && tfeb_id[1]!=11)) iKeepMe=0;	
	  if ((tfeb_id[1]==11 && tfeb_id[0]!=12) || (tfeb_id[1]==12 && tfeb_id[0]!=11)) iKeepMe=0;	
	}
	
	// apply hit threshold
	if (pestot<fHitThreshold) iKeepMe=0;
	else {
	// apply strip and sipm threshold
	std::vector<std::pair<int,float>> pes1 = tpesmap.find(tfeb_id[0])->second; 
	std::pair<int,float> ind_pes1=pes1[0];  // works only for simulation (the pes1 vector only has 2 elements)
	std::pair<int,float> ind_pes2=pes1[1];
	float tot1 = ind_pes1.second+ind_pes2.second;
	std::vector<std::pair<int,float>> pes2 = tpesmap.find(tfeb_id[1])->second; 
	std::pair<int,float> ind2_pes1=pes2[0];  // works only for simulation (the pes1 vector only has 2 elements)
	std::pair<int,float> ind2_pes2=pes2[1];
	float tot2 = ind2_pes1.second+ind2_pes2.second;
	if ( tot2<fStripThreshold || tot1<fStripThreshold ) iKeepMe=0;
	if (ind2_pes1.second < fSiPMThreshold || ind2_pes2.second<fSiPMThreshold) iKeepMe=0;
	if (ind_pes1.second < fSiPMThreshold || ind_pes2.second<fSiPMThreshold ) iKeepMe=0;
	//	if (iKeepMe==0) std::cout << "tot1 " << tot1 << " tot2 " << tot2 << std::endl;
	}
      } // if this is a MC hit
      if (iKeepMe) {
	
	// Create a corrected CRT hit
	crt::CRTHit crtHit = FillCrtHit(tfeb_id, tpesmap, pestot, time1,  time2,  time3,  time4,  time5, 
					plane, x, ex,y,ey,z,ez );
	
	CRTHitOutCol->push_back(crtHit);
	nHits++;
	if (fVerbose) std::cout << "hit created: time " << time5 << " x " <<  x << 
			" y " << y << " z " <<  z << std::endl;
      }  // keep this hit	    
    } // loop over hits
    
    event.put(std::move(CRTHitOutCol));

    if(fVerbose) std::cout<<"Number of CRT hits produced = "<<nHits<<std::endl;
      
    
  } // produce()
    
    void CRTSimHitCorr::endJob()
    {
      
    }
    
    crt::CRTHit CRTSimHitCorr::FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap, float peshit,double time1, double time2, double time3, double time4, double time5, 
int plane, double x, double ex, double y, double ey, double z, double ez){
	
	crt::CRTHit crtHit;
	crtHit.feb_id = tfeb_id;
	crtHit.pesmap = tpesmap;
	crtHit.peshit = peshit;
	crtHit.ts0_s = time1; 
	crtHit.ts0_s_corr = time2;
	crtHit.ts0_ns = time3;
	crtHit.ts0_ns_corr = time4;
	crtHit.ts1_ns = time5 ;
	crtHit.plane = plane;
	crtHit.x_pos = x;
	crtHit.x_err = ex;
	crtHit.y_pos = y; 
	crtHit.y_err = ey;
	crtHit.z_pos = z;
	crtHit.z_err = ez;
	return crtHit;
      }
      

      DEFINE_ART_MODULE(CRTSimHitCorr)

    }// namespace crt

namespace {


}
