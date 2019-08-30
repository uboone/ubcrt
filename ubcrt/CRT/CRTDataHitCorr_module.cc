////////////////////////////////////////////////////////////////////////////
/// Class:       CRTDataHitCorr
/// Module Type: producer
/// File:        CRTDataHitCorr_module.cc
///
/// Author:         Michelle Stancari
/// E-mail address: mstancar@fnal.gov
///
/////////////////////////////////////////////////////////////////////////////

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
  416.7, 644.7, 872.7, 1100.7, 188.7, 416.7, 644.7, 872.7, 1100.7, -230.0,                 //40-49
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

  
  class CRTDataHitCorr : public art::EDProducer {
  public:

    explicit CRTDataHitCorr(fhicl::ParameterSet const & p);

    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTDataHitCorr(CRTDataHitCorr const &) = delete;
    CRTDataHitCorr(CRTDataHitCorr &&) = delete;
    CRTDataHitCorr & operator = (CRTDataHitCorr const &) = delete; 
    CRTDataHitCorr & operator = (CRTDataHitCorr &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;

    void endJob() override;

    void reconfigure(fhicl::ParameterSet const & p);
    void CorrectAlignment(int plane, int feb1, int feb2, float &x, float &y, float &z);

    bool ApplyLightLossCorrection(int firstFEB,int secondFEB,float x, float y, float z, float &sf1, float &sf2);

    double ApplyT0Correction();

    crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t,std::vector<std::pair<int,float>>> tpesmap, 
			   float peshit, double time1, double time2, double time3, double time4, double time5, int plane,
			   double x, double ex, double y, double ey, double z, double ez); 




  private:

    // Params got from fcl file.......
    art::InputTag fCrtHitsIn_Label1;      ///< name of crt producer
    art::InputTag fCrtHitsIn_Label2;      ///< name of crt producer
    int fNumberCollections;
    bool fPlane3Only_Coll2;
    float fHitThreshold;
    float fStripThreshold;
    float fSiPMThreshold;

    //alignment params
    bool fRestorePE;
    bool fSumPE;
    bool fCorrectAlignment;
    bool fRemoveHits;
    bool fCorrectTiming;
    float fDistOffStrip;
			   //
    bool          fVerbose;             ///< print info
   
  }; // class CRTDataHitCorr
    
  CRTDataHitCorr::CRTDataHitCorr(fhicl::ParameterSet const & p)
  // Initialize member data here, if know don't want to reconfigure on the fly
  {
    // Call appropriate produces<>() functions here.
    produces< std::vector<crt::CRTHit> >();
    

    reconfigure(p);

  } // CRTDataHitCorr()

  void CRTDataHitCorr::reconfigure(fhicl::ParameterSet const & p)
  {
    fCrtHitsIn_Label1       = (p.get<art::InputTag> ("CrtHitsIn_Label1","merger")); 
    fCrtHitsIn_Label2       = (p.get<art::InputTag> ("CrtHitsIn_Label2","remerge")); 
    fNumberCollections      = (p.get<int>  ("NumberCollections",1));
    fPlane3Only_Coll2       = (p.get<bool> ("Plane3Only_Coll2",true));
    //
    fHitThreshold           = (p.get<float>("HitThreshold"  ,0.0));
    fStripThreshold         = (p.get<float>("StripThreshold",0.0));
    fSiPMThreshold          = (p.get<float>("SiPMThreshold" ,0.0));
    //alignment params
    fCorrectAlignment       = (p.get<bool> ("CorrectAlignment",true));
    fRestorePE              = (p.get<bool> ("RestorePE"       ,false));
    fSumPE                  = (p.get<bool> ("SumPE"           ,false));
    fRemoveHits             = (p.get<bool> ("RemoveHits"      ,false));
    fDistOffStrip           = (p.get<float>("DistOffStrip"    ,40.0));
    fVerbose                = (p.get<bool> ("Verbose"         ,false));
    fCorrectTiming          = (p.get<bool> ("CorrectTiming"   ,false));

  }

  void CRTDataHitCorr::beginJob()
    {
    if(fVerbose){std::cout<<"----------------- CRT Hit Reco Module -------------------"<<std::endl;}
    
  } // beginJob()


    double CRTDataHitCorr::ApplyT0Correction()
    {
      double correctionFactor = 0;
      //CRTDataQuality.... | crtt0Correction..... | ..................... | std::vector<anab::T0>....................... | ....1
      return correctionFactor;
      
    }

  void CRTDataHitCorr::CorrectAlignment(int plane, int feb1, int feb2, float &x, float &y, float &z)
  {
    if (plane==3) {x+=12.1;  y-=40.0;  z-=31.1;}
    else if (plane==2) {
      if  (feb1==15 ||  feb1==20 ||  feb1==46 ||  feb1==48 ||  feb1==50  ) y-=5.0;
      if  (feb2==15 ||  feb2==20 ||  feb2==46 ||  feb2==48 ||  feb2==50  ) y-=5.0;
      if  (feb1==16 ||  feb1==21 ||  feb1==47 ||  feb1==49 ||  feb1==51  ) y-=9.0;
      if  (feb2==16 ||  feb2==21 ||  feb2==47 ||  feb2==49 ||  feb2==51  ) y-=9.0;
      if  ((feb1>52 && feb1<56) || (feb2>52 && feb2<56)) y-=6.2;
      if ((feb1>31 && feb1<39) || (feb2>31 && feb2<39) ) z+=5.2;
      if ((feb1>38 && feb1<46) || (feb2>38 && feb2<46) ) z-=5.5;
    }
    else if (plane==1) {	  y+=20.7; z+=5.7;}
    else if (plane==0) {
      x-=4.4; 	  z+=32.8;
      if (feb1==11 || feb2==11) z-=10.4;
      if (feb1==12 || feb2==12) x-=4.9;
    }
  }
  

  bool CRTDataHitCorr::ApplyLightLossCorrection(int firstFEB,int secondFEB,float x, float y, float z, float &sf1, float &sf2)
  {
    // If requested, calculate corrections to hit pe values for light loss to attenuation 
    // If requested,  remove hits whose 3D psoition is far from the actual scintillator
    //     The latter should be done before the 3D position is adjusted for alignment??
    
    float b=1085.0; // cm
    int m1 = feb2mod[firstFEB]; int m2=feb2mod[secondFEB];
    float len1=mod2length[m1]; float len2=mod2length[m2]; 
    int or1=mod2orient[m1]; int or2=mod2orient[m2];
    float sipm_pos1=sipm_pos[m1]; float sipm_pos2=sipm_pos[m2];
    int sipm_end1 = mod2end[m1];  int sipm_end2 = mod2end[m2];
    float disttoreadout1 = 0.0;
    if (or1==0) disttoreadout1=sipm_end1*(x-sipm_pos1);
    else if (or1==1)  disttoreadout1=sipm_end1*(y-sipm_pos1);
    else  disttoreadout1=sipm_end1*(z-sipm_pos1);
    //
    if (fRemoveHits) {
      if (disttoreadout1<-1.0*fDistOffStrip) return false;
      else if (disttoreadout1>(len1+fDistOffStrip)) return false;
    }
    if (disttoreadout1>len1) disttoreadout1=len1;
    else if (disttoreadout1<0) disttoreadout1=0.0;
    if (fRestorePE) { // correct hit pe for light losses to attenuation
      sf1=pow(disttoreadout1+b,2.0)/b/b;
    }
    float disttoreadout2 = 0.0;
    if (or2==0) disttoreadout2=sipm_end2*(x-sipm_pos2);
    else if (or2==1)  disttoreadout2=sipm_end2*(y-sipm_pos2);
    else  disttoreadout2=sipm_end2*(z-sipm_pos2);
    if (fRemoveHits) {
      if (disttoreadout2<-1.0*fDistOffStrip) return false;
      else if (disttoreadout2>(len2+fDistOffStrip)) return false;
    }
    if (disttoreadout2>len2) disttoreadout2=len2;
    else if (disttoreadout2<0) disttoreadout2=0.0;
    if (fRestorePE) { // correct hit pe for light losses to attenuation
      sf2=pow(disttoreadout2+b,2.0)/b/b;
    }	  // if RestorePE is turned on 
    return true;
  }
  
    
  void CRTDataHitCorr::produce(art::Event & event)
  {
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    // Place to store corrected CRThits as they are created
    std::unique_ptr<std::vector<crt::CRTHit>> CRTHitOutCol( new std::vector<crt::CRTHit>);


    std::vector<crt::CRTHit> crtHitInList;
    // Retrieve first list of CRT hits
    art::Handle< std::vector<crt::CRTHit>> crtHitsInHandle;
    event.getByLabel(fCrtHitsIn_Label1, crtHitsInHandle);
    //check to make sure the data we asked for is valid
    if(!crtHitsInHandle.isValid()){
      std::cout << "Run " << event.run() << ", subrun " << event.subRun()
		<< ", event " << event.event() << " has zero"
		<< " CRTHits " << " in module " << fCrtHitsIn_Label1 << std::endl;
      std::cout << std::endl;
      //add protection here
      event.put(std::move(CRTHitOutCol));
      return;
    }

    std::vector<crt::CRTHit> const& crtHitInList1(*crtHitsInHandle);
    crtHitInList.insert(crtHitInList.end(), crtHitInList1.begin(), crtHitInList1.end());
    if(fVerbose) std::cout<<"Number of CRT hits read in= "<<crtHitInList.size()<< 
		   " after first collection" << std::endl;
    //
    if (fNumberCollections>1) {
    // Retrieve second list of CRT hits
    event.getByLabel(fCrtHitsIn_Label2, crtHitsInHandle);
    //check to make sure the data we asked for is valid
    if(!crtHitsInHandle.isValid()){
      std::cout << "Run " << event.run() << ", subrun " << event.subRun()
		<< ", event " << event.event() << " has zero"
		<< " CRTHits " << " in module " << fCrtHitsIn_Label2 << std::endl;
      std::cout << " Skipping this CRT hit collection " << std::endl;
    }
    std::vector<crt::CRTHit> const& crtHitInList2(*crtHitsInHandle);
    if (fPlane3Only_Coll2) {
      for (size_t i = 0; i < crtHitInList2.size(); i++){
	if (crtHitInList2[i].plane==3)  crtHitInList.insert(crtHitInList.end(), crtHitInList2[i]);
      }
    }
    else     crtHitInList.insert(crtHitInList.end(), crtHitInList2.begin(), crtHitInList2.end());
    if(fVerbose) std::cout<<"Number of CRT hits read in= "<<crtHitInList.size()<< " after second collection" << std::endl;
    }

    if(fVerbose) std::cout<<"In data correction "<<crtHitInList.size()<<" "<<fCrtHitsIn_Label1<<"\n";

    for (size_t i = 0; i < crtHitInList.size(); i++){
      // Let's start positive, let's keep this hit
      bool keepMe = true;

      //Get hit
      crt::CRTHit thisCrtHit = crtHitInList[i];

      
      // Get all memebers you need to check if this is data or MC.
      // only change/remove MC hits. 
      std::vector<uint8_t> tfeb_id = thisCrtHit.feb_id;
      int firstFEB  = (int)tfeb_id[0];
      int secondFEB = (int)tfeb_id[1];
      std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap=thisCrtHit.pesmap;
      std::vector<std::pair<int,float>> firstStrip  = tpesmap.find(firstFEB)->second; 
      std::vector<std::pair<int,float>> secondStrip = tpesmap.find(secondFEB)->second; 
      // if this is mc, push the hit and continue with the next
      if (firstStrip.size() == 2) { CRTHitOutCol->push_back(thisCrtHit); continue;}
      // At this point, throw an exeption if the size is strange
      if (firstStrip.size() !=  32) {std::cerr << "\033[93m[ERROR]\033[00m Hit has wrong number of strips: "<<tpesmap.size()  << std::endl;
	event.put(std::move(CRTHitOutCol));
	return; }
      
      // Ok, now we just made sure we have only data hit... let modify them!
      float time1 = thisCrtHit.ts0_s;
      float time2 = thisCrtHit.ts0_s_corr;
      float time3 = thisCrtHit.ts0_ns;
      float time4 = thisCrtHit.ts0_ns_corr;
      float time5 = thisCrtHit.ts1_ns;
      int   plane = thisCrtHit.plane;
      float    x  = thisCrtHit.x_pos;
      float    ex = thisCrtHit.x_err;
      float    y  = thisCrtHit.y_pos;
      float    ey = thisCrtHit.y_err;
      float    z  = thisCrtHit.z_pos;
      float    ez = thisCrtHit.z_err;
      float pestot = thisCrtHit.peshit;      
      

      // remove if under threshold
      if (pestot<fHitThreshold) continue;

      //Calculate attenuation factor 
      float sf1=1.0; float sf2=1.0;	
      if (fRemoveHits || fRestorePE)  keepMe = ApplyLightLossCorrection(firstFEB,secondFEB,x,y,z,sf1,sf1);
      if (!keepMe) continue;

      // apply strip and sipm threshold
      std::pair<int,float> ind_pes1,ind_pes2,ind2_pes1,ind2_pes2;
      float pmax1=0.0; float pmax2=0.0;
      float totA=0; float totB=0;
      std::vector<std::pair<int,float>> pesAnew;  pesAnew.clear();
      std::vector<std::pair<int,float>> pesBnew;  pesBnew.clear();
      for (int ii=0;ii<32;ii+=2){ 
	std::pair<int,float> ind_pesa=firstStrip[ii];  
	std::pair<int,float> ind_pesb=firstStrip[ii+1];
	// find strip with largest pe value for applying threshold
	if ((ind_pesa.second+ind_pesb.second)>pmax1) {
	  ind_pes1=ind_pesa; ind_pes2=ind_pesb;
	  pmax1=ind_pesa.second+ind_pesb.second;	    
	}
	totA+=ind_pesa.second+ind_pesb.second;	    
	std::pair<int,float> pesA0(ind_pesa.first,sf1*ind_pesa.second);
	std::pair<int,float> pesA1(ind_pesb.first,sf1*ind_pesb.second);
	pesAnew.push_back(pesA0); pesAnew.push_back(pesA1);
	//
	ind_pesa=secondStrip[ii];  
	ind_pesb=secondStrip[ii+1];
	if ((ind_pesa.second+ind_pesb.second)>pmax2) {
	  ind2_pes1=ind_pesa; ind2_pes2=ind_pesb;
	  pmax2=ind_pesa.second+ind_pesb.second;	    
	}
	totB+=ind_pesa.second+ind_pesb.second;	    
	std::pair<int,float> pesB0(ind_pesa.first,sf2*ind_pesa.second);
	std::pair<int,float> pesB1(ind_pesb.first,sf2*ind_pesb.second);
	pesBnew.push_back(pesB0); pesBnew.push_back(pesB1);
      }

	
	// issue here reconstructing exactly the report peshit for the hit - hits are created using only the 
	//    strip with the max ADC.   But strip values are saved in units of pe, not ADC.  We take the strip
	//    with the max pe
	if ( pmax2<fStripThreshold || pmax2<fStripThreshold ) continue;
	if (ind2_pes1.second < fSiPMThreshold || ind2_pes2.second<fSiPMThreshold) continue;
	if (ind_pes1.second < fSiPMThreshold || ind_pes2.second<fSiPMThreshold ) continue;
	//
	if (fSumPE) pestot=totA*sf1+totB*sf2;
	else pestot=pmax1*sf1+pmax2*sf2;

	tpesmap.clear();
	tpesmap.emplace(firstFEB,pesAnew);
	tpesmap.emplace(secondFEB,pesBnew);
	

	if (fCorrectAlignment)  CorrectAlignment(plane, firstFEB, secondFEB, x, y, z);
	  
	if (fCorrectTiming ) {
	  double t0correction = ApplyT0Correction();
	  time3 += t0correction; //<-- wrong
	  time4 += t0correction; //<-- wrong
	}

	// Create a corrected CRT hit
	crt::CRTHit crtHit = FillCrtHit(tfeb_id, tpesmap, pestot, time1,  time2,  time3,  time4,  time5, 
					plane, x, ex,y,ey,z,ez );
	
	CRTHitOutCol->push_back(crtHit);
	
    } // loop over hits
    
    event.put(std::move(CRTHitOutCol));

      
    
  } // produce()
    
    void CRTDataHitCorr::endJob()
    {
      
    }
    
    crt::CRTHit CRTDataHitCorr::FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap, float peshit,double time1, double time2, double time3, double time4, double time5, 
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
      

      DEFINE_ART_MODULE(CRTDataHitCorr)

    }// namespace crt

namespace {


}
