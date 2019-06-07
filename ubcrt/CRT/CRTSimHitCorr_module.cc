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
#include <fstream>
#include <map>
#include <iterator>
#include <algorithm>
#include <vector>
#include <unordered_set>


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
    
   
    void DBCall( std::vector<std::pair<int,int>> &deadMap );
    bool isHitFromDeadChannels(int febNumber1, int channel1Number1, int channel1Number2 , 
			       int febNumber2, int channel2Number1, int channel2Number2 , std::vector<std::pair<int,int>> deadMap );
 
    crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t,std::vector<std::pair<int,float>>> tpesmap, 
			   float peshit, double time1, double time2, double time3, double time4, double time5, int plane,
			   double x, double ex, double y, double ey, double z, double ez); 
    
    void ScaleMCtime  (float &time1,float &time5) {time1/=5.;time5/=5.;};
    

    bool BottomSections(int firstFEB, int secondFEB);
    bool PipeSections(int firstFEB, int secondFEB);			   
    bool FTSections(int firstFEB, int secondFEB);			   
    bool TopSections(int firstFEB, int secondFEB);	
    bool ApplyDetectorResponse(crt::CRTHit thisCrtHit,  std::map<uint8_t, std::vector<std::pair<int,float>>> &tpesmap, float &pestot );
    /*

    void RestorePE             (crt::CRTHit hit);
    void ApplyDetectorResponse (crt::CRTHit hit);
    void RemoveHits            (crt::CRTHit hit);
    void MaskDeadChannels      (crt::CRTHit hit);
    void SimulatedSaturation   (crt::CRTHit hit);
			   */

    
  private:
    // Params got from fcl file.......
    art::InputTag fCrtHitsIn_Label;     ///< name of crt producer
    bool  fScaleMCtime;                 ///< turns off/on bug fix for hit times
    float fHitThreshold;
    float fStripThreshold;
    float fSiPMThreshold;
    float fPEscaleFactor;
    float fElectNoise;
    float fDistOffStrip;
    bool  fRestorePE;
    bool  fRemoveBottomHits;
    bool  fApplyDetectorResponse;
    bool  fVerbose;
    bool  fRemoveHits;
    bool  fMaskDeadChannels;
    std::vector<int> fDeadFEB; 
    std::vector<int> fDeadChannels;
    bool  fSections;
    bool  fSimulatedSaturation;
    
    ///< print info
    CLHEP::HepRandomEngine& fEngine;

    std::vector<int> MichelleFEB      {29, 30, 32, 37, 37, 37, 37, 38, 41, 46, 109, 111, 113, 113,  124};   //117,
    std::vector<int> MichelleChannels {23,  1,  3,  2,  6, 12, 26,  8,  0,  7,  11,  31,   6,   8,   14};   // 21, 
			   
  }; // class CRTSimHitCorr
    

    CRTSimHitCorr::CRTSimHitCorr(fhicl::ParameterSet const & p) : EDProducer{p},     
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

      // Default parameters
      fCrtHitsIn_Label       = (p.get<art::InputTag> ("CrtHitsIn_Label"      ,"crtsimhit")); 
      fScaleMCtime           = (p.get<bool>          ("ScaleMCtime"          ,false));
      fHitThreshold          = (p.get<float>         ("HitThreshold"         ,0.0));
      fStripThreshold        = (p.get<float>         ("StripThreshold"       ,7.75));
      fSiPMThreshold         = (p.get<float>         ("SiPMThreshold"        ,0.0));
      fPEscaleFactor         = (p.get<float>         ("PEscaleFactor"        ,1.525));
      fElectNoise            = (p.get<float>         ("ElectNoise"           ,2.0));
      fRestorePE             = (p.get<bool>          ("RestorePE"            ,true));
      fRemoveHits            = (p.get<bool>          ("RemoveHits"           ,true));
      fRemoveBottomHits      = (p.get<bool>          ("RemoveBottomHits"     ,true));
      fApplyDetectorResponse = (p.get<bool>          ("ApplyDetectorResponse",true));
      fMaskDeadChannels      = (p.get<bool>          ("MaskDeadChannels"     ,true));
      fSections              = (p.get<bool>          ("Sections"             ,true));
      fSimulatedSaturation   = (p.get<bool>          ("SimulatedSaturation"  ,true));

      fDeadFEB               = (p.get< std::vector<int> > ("DeadFEB"    ,  MichelleFEB ));
      fDeadChannels          = (p.get< std::vector<int> > ("DeadChannel",  MichelleChannels ));

      fVerbose               = (p.get<bool>          ("Verbose"              ,false));
    }

    void CRTSimHitCorr::beginJob()
    {
      if(fVerbose){std::cout<<"----------------- CRT Hit Correction Module -------------------"<<std::endl;}
    } // beginJob()


    void CRTSimHitCorr::endJob()
    {
    }

    void CRTSimHitCorr::DBCall( std::vector<std::pair<int,int>> &deadMap )
    {
      deadMap.clear();
      if (fDeadFEB.size() != fDeadChannels.size() ) {
	std::cout<<"FEB and Channel sizes are different. I am not applying the dead channel masking\n";
	return;
      }
      
      for (size_t i = 0; i < fDeadFEB.size(); i++ ){
	std::pair<int,int> p1 (fDeadFEB[i], fDeadChannels[i]);
	deadMap.push_back(p1);
      }
    }

    bool  CRTSimHitCorr::isHitFromDeadChannels(int febNumber1, int channel1Number1, int channel1Number2 , 
					       int febNumber2, int channel2Number1, int channel2Number2, std::vector<std::pair<int,int>> deadMap )
    {
      for(auto const& value: deadMap)
	{ 
	  if (value.first == febNumber1 && value.second == channel1Number1) {return true;}
	  if (value.first == febNumber1 && value.second == channel1Number2) {return true;}
	  if (value.first == febNumber2 && value.second == channel2Number1) {return true;}
	  if (value.first == febNumber2 && value.second == channel2Number2) {return true;}
	}
      return false;
    }

  crt::CRTHit CRTSimHitCorr::FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap, float peshit,
					double time1, double time2, double time3, double time4, double time5, 
					int plane, double x, double ex, double y, double ey, double z, double ez) {
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

  bool CRTSimHitCorr::BottomSections(int firstFEB, int secondFEB)
  {    
    std::unordered_set<int> sectionA = {11,12};
    std::unordered_set<int> sectionB = {14,17,23,24};
    std::unordered_set<int> sectionC = {18,19,22,23};

    // Check if both FEBs are in section A
    const bool first_inA  = sectionA.find(firstFEB)  != sectionA.end();
    const bool second_inA = sectionA.find(secondFEB) != sectionA.end();
    // Check if both FEBs are in section B
    const bool first_inB  = sectionB.find(firstFEB)  != sectionB.end();
    const bool second_inB = sectionB.find(secondFEB) != sectionB.end();

    // Check if both FEBs are in section C
    const bool first_inC  = sectionC.find(firstFEB)  != sectionC.end();
    const bool second_inC = sectionC.find(secondFEB) != sectionC.end();
	  
    if ( !((first_inA && second_inA) ||  (first_inB && second_inB) ||  (first_inC && second_inC) ) ) return false;
    return true;
  }


  bool CRTSimHitCorr::FTSections(int firstFEB, int secondFEB)
  {    
    std::unordered_set<int> sectionA = {52,31,29,60,61};
    std::unordered_set<int> sectionB = {29,28,27,58,59};
    std::unordered_set<int> sectionC = {27,26,30,56,57};

    // Check if both FEBs are in section A
    const bool first_inA  = sectionA.find(firstFEB)  != sectionA.end();
    const bool second_inA = sectionA.find(secondFEB) != sectionA.end();
    // Check if both FEBs are in section B
    const bool first_inB  = sectionB.find(firstFEB)  != sectionB.end();
    const bool second_inB = sectionB.find(secondFEB) != sectionB.end();

    // Check if both FEBs are in section C
    const bool first_inC  = sectionC.find(firstFEB)  != sectionC.end();
    const bool second_inC = sectionC.find(secondFEB) != sectionC.end();
	  
    if ( !((first_inA && second_inA) ||  (first_inB && second_inB) ||  (first_inC && second_inC) ) ) return false;
    return true;
  }


  bool CRTSimHitCorr::TopSections(int firstFEB, int secondFEB)
  {    
    std::unordered_set<int> sectionA = {109,105,195,123,113,114,115};
    std::unordered_set<int> sectionB = {106,124,107,108,116,117,118,127};
    std::unordered_set<int> sectionC = {128,125,126,111,112,119,120,121,129};

    // Check if both FEBs are in section A
    const bool first_inA  = sectionA.find(firstFEB)  != sectionA.end();
    const bool second_inA = sectionA.find(secondFEB) != sectionA.end();
    // Check if both FEBs are in section B
    const bool first_inB  = sectionB.find(firstFEB)  != sectionB.end();
    const bool second_inB = sectionB.find(secondFEB) != sectionB.end();

    // Check if both FEBs are in section C
    const bool first_inC  = sectionC.find(firstFEB)  != sectionC.end();
    const bool second_inC = sectionC.find(secondFEB) != sectionC.end();
	  
    if ( !((first_inA && second_inA) ||  (first_inB && second_inB) ||  (first_inC && second_inC) ) ) return false;
    return true;
  }
					
  bool CRTSimHitCorr::PipeSections(int firstFEB, int secondFEB)
  {
    std::unordered_set<int> sectionA = {53,37,33,34,15,16,20,21};
    std::unordered_set<int> sectionB = {54,34,35,36,20,21,46,47,48,49};
    std::unordered_set<int> sectionC = {55,36,32,38,46,47,48,49,50,51};

    std::unordered_set<int> sectionD = {39,40,41,15,16,20,21,46,47};
    std::unordered_set<int> sectionE = {41,42,43,46,47,48,49};    
    std::unordered_set<int> sectionF = {43,44,45,48,49,50,51};

    //if (firstFEB==40||secondFEB==40) std::cout<<secondFEB<<" "<<firstFEB<<"\n";
    // Check if both FEBs are in section A
    const bool first_inA  = sectionA.find(firstFEB)  != sectionA.end();
    const bool second_inA = sectionA.find(secondFEB) != sectionA.end();
    // Check if both FEBs are in section B
    const bool first_inB  = sectionB.find(firstFEB)  != sectionB.end();
    const bool second_inB = sectionB.find(secondFEB) != sectionB.end();
    // Check if both FEBs are in section C
    const bool first_inC  = sectionC.find(firstFEB)  != sectionC.end();
    const bool second_inC = sectionC.find(secondFEB) != sectionC.end();
    // Check if both FEBs are in section D
    const bool first_inD  = sectionD.find(firstFEB)  != sectionD.end();
    const bool second_inD = sectionD.find(secondFEB) != sectionD.end();
    // Check if both FEBs are in section E
    const bool first_inE  = sectionE.find(firstFEB)  != sectionE.end();
    const bool second_inE = sectionE.find(secondFEB) != sectionE.end();
    // Check if both FEBs are in section F
    const bool first_inF  = sectionF.find(firstFEB)  != sectionF.end();
    const bool second_inF = sectionF.find(secondFEB) != sectionF.end();
    /*
    if ( !first_inA &&  !first_inB &&  !first_inC &&
	 !first_inD &&  !first_inE &&  !first_inF     ) {std::cout<<firstFEB<<"\n";}

    if ( !second_inA &&  !second_inB &&  !second_inC &&
	 !second_inD &&  !second_inE &&  !second_inF     ) {std::cout<<secondFEB<<"\n";}
    */
    if ( !((first_inA && second_inA) ||  (first_inB && second_inB) ||  (first_inC && second_inC) ||
	   (first_inD && second_inD) ||  (first_inE && second_inE) ||  (first_inF && second_inF)    ) ) {return false;}
    return true;
  }

    
    void CRTSimHitCorr::produce(art::Event & event)
    {
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
      // Load the crt dead channel map only once per event
      std::vector<std::pair<int,int>> deadMap;
      if (fMaskDeadChannels){
	DBCall(deadMap);
      }

      std::vector<crt::CRTHit> const& crtHitInList(*crtHitsInHandle);
      if(fVerbose) std::cout<<"Number of CRT hits read in= "<<crtHitInList.size()<< std::endl;
      

      //Loop on CRT Hits
      for (size_t i = 0; i < crtHitInList.size(); i++){
	// Let's start positive, let's keep this hit
	bool keepMe = 1;

	// Get hit
	crt::CRTHit thisCrtHit = crtHitInList[i];
	// Get all memebers you need to check if this is data or MC.
	// only change/remove MC hits. 
	std::vector<uint8_t> tfeb_id = thisCrtHit.feb_id;
	int firstFEB  = (int)tfeb_id[0];
	int secondFEB = (int)tfeb_id[1];
	std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap=thisCrtHit.pesmap;
	std::vector<std::pair<int,float>> firstStrip  = tpesmap.find(firstFEB)->second; 
	std::vector<std::pair<int,float>> secondStrip = tpesmap.find(secondFEB)->second; 
	// if this is data, push the hit and continue with the next
	if (firstStrip.size() == 32) { CRTHitOutCol->push_back(thisCrtHit); continue;}
	// At this point, throw an exeption if the size is strange
	if (firstStrip.size() !=  2) {std::cerr << "\033[93m[ERROR]\033[00m Hit has wrong number of strips: "<<tpesmap.size()  << std::endl;
	  event.put(std::move(CRTHitOutCol));
	  return; }
	
	// Ok, now we just made sure we have only MC hit... let modify them!
	
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


	// Modify the timing of ts1 due to bug in
	// in Oct2018 simulation. May not apply anymore
	if (fScaleMCtime) ScaleMCtime(time1,time5);


	// Would you like restrict the strip crossing to overlapping modules?
	if (fSections)
	  {
	    if (plane == 0) keepMe = BottomSections(firstFEB, secondFEB);
	    else if (plane == 1) keepMe = FTSections(firstFEB, secondFEB);
	    else if (plane == 2) keepMe = PipeSections(firstFEB, secondFEB);
	    else if (plane == 3) keepMe = TopSections(firstFEB, secondFEB);
	  }
	if (!keepMe) continue; // If this hit is to trash, skip the rest
	

	//Let's mask the dead channels
	if (fMaskDeadChannels) keepMe = !(isHitFromDeadChannels(firstFEB , firstStrip [0].first, firstStrip [0].first, 
								secondFEB, secondStrip[0].first, secondStrip[0].first, 
								deadMap));

	if (!keepMe) continue; // If this hit is to trash, skip the rest

	// Apply the detector response & threshold at the sipm & strip level
	if (fApplyDetectorResponse) { keepMe = ApplyDetectorResponse(thisCrtHit,  tpesmap, pestot); }
	// after that, check the hit threshold
	if (pestot<fHitThreshold) continue;

	
	if (keepMe) {
	  // Create a corrected CRT hit
	  crt::CRTHit crtHit = FillCrtHit(tfeb_id, tpesmap, pestot,
					  time1,  time2,  time3,  time4,  time5, 
					  plane, x, ex,y,ey,z,ez );
	  CRTHitOutCol->push_back(crtHit);
	}// If the hit is to be keptSW
      }// End of crt hits loop
      event.put(std::move(CRTHitOutCol));
    }// End of the produce function



  bool CRTSimHitCorr::ApplyDetectorResponse(crt::CRTHit thisCrtHit,  std::map<uint8_t , std::vector<std::pair<int,float>>>  &tpesmap, float &pestot )
  {
    // Get the strips info!
    std::vector<uint8_t> tfeb_id = thisCrtHit.feb_id;
    int first_feb  = (int)tfeb_id[0];
    
    if (fPEscaleFactor<0) fPEscaleFactor=1.0; 
    
    // ---------------------- first strip -------------------------------
    double distToReadout=0;
    int this_mod = feb2mod[first_feb];
    if (this_mod<0 || this_mod>72) { std::cout << "bad module number for feb " << first_feb << std::endl; return false;}
    // Calculate the distance between hit and sipm... it depends on the orientation of the module (= 3 possibility)
    if      (mod2orient[this_mod]==0) distToReadout=mod2end[this_mod]*(thisCrtHit.x_pos-sipm_pos[this_mod]);
    else if (mod2orient[this_mod]==1) distToReadout=mod2end[this_mod]*(thisCrtHit.y_pos-sipm_pos[this_mod]);
    else                              distToReadout=mod2end[this_mod]*(thisCrtHit.z_pos-sipm_pos[this_mod]);
    
    // And fudge a bit with the ends
    if (fRemoveHits) {
      if (distToReadout<-1.0*fDistOffStrip) return false;
      else if (distToReadout>(mod2length[this_mod]+fDistOffStrip)) return false;
    }
    if (distToReadout>mod2length[this_mod]) {
      distToReadout=mod2length[this_mod];
    }
    else if (distToReadout<0) distToReadout=0.0;
    
    
    // Setup the strip width (the top is wider)
    double stripWidth = 10.8;
    if (thisCrtHit.plane == 3) stripWidth = 11.2;

    // This is where the fun starts
    double b = 1085.0;
    double pe_sf_A = b*b/pow(distToReadout+b,2.0);
    
    // Take the crt hit map, 
    // Find the first feb for this hit
    std::map<uint8_t, std::vector<std::pair<int,float>>> thispesmap=thisCrtHit.pesmap;
    std::vector<std::pair<int,float>> pesA  = thispesmap.find(first_feb)->second; 
    // Take the signal on the single sipm and scale it by 2 factors: fPEscaleFactor and pe_sf_A
    float pesA_0=(fPEscaleFactor*pe_sf_A)*pesA[0].second;
    float pesB_0=(fPEscaleFactor*pe_sf_A)*pesA[1].second;
    // Calculate the total pes
    float pestot_0=pesA_0+pesB_0;
    // Calculate the distance in the "across" direction 
    float distA=(pesA_0/pestot_0)*stripWidth;  // 10.8 is strip width in cm for side, 11.2 for the top
    float distB=(pesB_0/pestot_0)*stripWidth;  // 10.8 is strip width in cm for side, 11.2 for the top
    // Calculate attenuation in the "across" direction 
    float absA=exp(-1.0*distA/8.5); // 
    float absB=exp(-1.0*distB/8.5); //
    // Ridistribute the signal on single sipm considering the attenuation
    pesA_0=pestot_0*(absA)/(absA+absB);
    pesB_0=pestot_0*(absB)/(absA+absB);
    // smear these values to simulate the number of produced photoelectrons
    float npe0 = CLHEP::RandPoisson::shoot(&fEngine, pesA_0);
    float npe1 = CLHEP::RandPoisson::shoot(&fEngine, pesB_0);
    // SiPM and ADC response: Npe to ADC counts
    float pesA_sm = CLHEP::RandGauss::shoot(&fEngine, npe0, fElectNoise * sqrt(npe0)); // fElectNoise guess 0.085
    float pesB_sm = CLHEP::RandGauss::shoot(&fEngine, npe1, fElectNoise * sqrt(npe0)); // fElectNoise guess 0.085
    if (fVerbose)std::cout<<"DEBUGG  "<< fElectNoise<<"  "<<pesA_0<<"   "<<npe0<<" --   "<<pesA_sm<<"  ";
	  
    if (fSimulatedSaturation)
      {
	// If the pe is too high, the sipm saturated and returns only the maximum pe.
	// We estimate this saturation from data in the following way:
	// For a 12 bit adc you have a range of 4095, the average gain in data (i.e. ADC/pe conversion) is 
	// 40 ADC/pe, so 4095/40 = 102.375 is the max pe recorded... let's put a cap on that!
	if (pesA_sm > 102.375 ) pesA_sm = 102.375;
	if (pesB_sm > 102.375 ) pesB_sm = 102.375;
      }

    if (fRestorePE) {
      float sf2=pow(distToReadout+b,2.0)/b/b;
      pesA_sm *= sf2;
      pesB_sm *= sf2;
    }
    
    if (pesA_sm< fSiPMThreshold) return false;
    if (pesB_sm< fSiPMThreshold) return false;
    if ( (pesA_sm+pesB_sm) <fStripThreshold) return false;

    // Put these values back into the single simulated sipms
    std::pair<int,float> pesA0(pesA[0].first,pesA_sm);
    std::pair<int,float> pesA1(pesA[1].first,pesB_sm);
    std::vector<std::pair<int,float>> pesAnew;
    pesAnew.push_back(pesA0); pesAnew.push_back(pesA1);

    tpesmap.clear();
    tpesmap.emplace(tfeb_id[0],pesAnew);
    // Calculate the contribution of the first stript to the total pe of the hit
    pestot = pesA_sm+pesB_sm;


    // ---------------------- second strip -------------------------------
    int second_feb = (int)tfeb_id[1];
    std::vector<std::pair<int,float>> pesB  = thispesmap.find(second_feb)->second; 
    distToReadout=0;
    this_mod = feb2mod[second_feb];
    if (this_mod<0 || this_mod>72) 
      std::cout << "bad module number for feb " << second_feb << std::endl;
    else {	    
      if (mod2orient[this_mod]==0) distToReadout=mod2end[this_mod]*(thisCrtHit.x_pos-sipm_pos[this_mod]);
      else if (mod2orient[this_mod]==1) distToReadout=mod2end[this_mod]*(thisCrtHit.y_pos-sipm_pos[this_mod]);
      else distToReadout=mod2end[this_mod]*(thisCrtHit.z_pos-sipm_pos[this_mod]);
	    
      if (fRemoveHits) {
	if (distToReadout<-1.0*fDistOffStrip) return false;
	else if (distToReadout>(mod2length[this_mod]+fDistOffStrip)) return false;
      }
      if (distToReadout>mod2length[this_mod]) {
	distToReadout=mod2length[this_mod];
      }
      else if (distToReadout<0) distToReadout=0.0;
    }	  
    double pe_sf_B = b*b/pow(distToReadout+b,2.0);
    pesA_0=pe_sf_B*fPEscaleFactor*pesB[0].second;
    pesB_0=pe_sf_B*fPEscaleFactor*pesB[1].second;
    pestot_0=pesA_0+pesB_0;
    distA=(pesA_0/pestot_0)*stripWidth;  // 10.8 is strip width in cm for side, 11.2 for the top
    distB=(pesB_0/pestot_0)*stripWidth;  // 10.8 is strip width in cm for side, 11.2 for the top
    absA=exp(-1.0*distA/8.5);
    absB=exp(-1.0*distB/8.5);
    pesA_0=pestot_0*(absA)/(absA+absB);
    pesB_0=pestot_0*(absB)/(absA+absB);
    // smear these values
    npe0 = CLHEP::RandPoisson::shoot(&fEngine, pesA_0);
    npe1 = CLHEP::RandPoisson::shoot(&fEngine, pesB_0);
    // SiPM and ADC response: Npe to ADC counts
    pesA_sm = CLHEP::RandGauss::shoot(&fEngine, npe0, fElectNoise * sqrt(npe0));
    pesB_sm = CLHEP::RandGauss::shoot(&fEngine, npe1, fElectNoise * sqrt(npe1));
    if (fSimulatedSaturation)
      {
	// If the pe is too high, the sipm saturated and returns only the maximum pe.
	// We estimate this saturation from data in the following way:
	// For a 12 bit adc you have a range of 4095, the average gain in data (i.e. ADC/pe conversion) is 
	// 40 ADC/pe, so 4095/40 = 102.375 is the max pe recorded... let's put a cap on that!
	if (pesA_sm > 102.375 ) pesA_sm = 102.375;
	if (pesB_sm > 102.375 ) pesB_sm = 102.375;
      }
    if (fRestorePE) {
      float sf2=pow(distToReadout+b,2.0)/b/b;
      pesA_sm *= sf2;
      pesB_sm *= sf2;
    }
    if (pesA_sm< fSiPMThreshold) return false;
    if (pesB_sm< fSiPMThreshold) return false;
    if ( (pesA_sm+pesB_sm) <fStripThreshold) return false;

    std::pair<int,float> pesB0(pesB[0].first,pesA_sm);
    std::pair<int,float> pesB1(pesB[1].first,pesB_sm);
    std::vector<std::pair<int,float>> pesBnew;
    pesBnew.push_back(pesB0); pesBnew.push_back(pesB1);
    tpesmap.emplace(tfeb_id[1],pesBnew);
    pestot += pesA_sm+pesB_sm;

    return true;
  }

	
					    	


  
  DEFINE_ART_MODULE(CRTSimHitCorr)
}// namespace crt


     
