////////////////////////////////////////////////////////////////////////////
/// Class:        CRTOverlayHotFix
/// Module Type: producer
/// File:         CRTOverlayHotFix_module.cc
///
/// Author:         Elena Gramellini
/// E-mail address: elenag@fnal.gov
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

namespace crt{
  
  class  CRTOverlayHotFix : public art::EDProducer {
  public:
    explicit  CRTOverlayHotFix(fhicl::ParameterSet const & p);
    
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    CRTOverlayHotFix( CRTOverlayHotFix const &) = delete;
    CRTOverlayHotFix( CRTOverlayHotFix &&) = delete;
    CRTOverlayHotFix & operator = ( CRTOverlayHotFix const &) = delete; 
    CRTOverlayHotFix & operator = ( CRTOverlayHotFix &&) = delete;
    
    // Required functions.
    void produce(art::Event & e) override;
    
    // Selected optional functions.
    void beginJob() override;
    void endJob() override;
    void reconfigure(fhicl::ParameterSet const & p);
    
    crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap, 
			   float peshit, double time1, double time2, double time3, double time4, double time5, int plane,
			   double x, double ex, double y, double ey, double z, double ez); 
    
  private:
    // Params got from fcl file.......
    art::InputTag fCrtHitsIn_Label;      ///< name of crt producer
    bool  fVerbose;             ///< print info
    bool  fCheckHitGeo;             ///< print info
  }; // class  CRTOverlayHotFix
  
  
  CRTOverlayHotFix:: CRTOverlayHotFix(fhicl::ParameterSet const & p)
  // Initialize member data here, if know don't want to reconfigure on the fly
  {
    // Call appropriate produces<>() functions here.
    produces< std::vector<crt::CRTHit> >();
    reconfigure(p);
  } //  CRTOverlayHotFix()
  
  
  void  CRTOverlayHotFix::reconfigure(fhicl::ParameterSet const & p)
  {
    fCrtHitsIn_Label   = (p.get<art::InputTag> ("CrtHitsIn_Label","mixer")); 
    fVerbose            = (p.get<bool> ("Verbose",true));
  }

  void  CRTOverlayHotFix::beginJob()
  {
  } // beginJob()
    
  void  CRTOverlayHotFix::produce(art::Event & event)
  {
    int nHits = 0;
    
    if(fVerbose)  std::cout<<"============================================"<<std::endl
			   <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
			   <<"============================================"<<std::endl;

    
    // Place to store corrected CRThits as they are created
    std::unique_ptr<std::vector<crt::CRTHit>> CRTHitOutCol( new std::vector<crt::CRTHit>);
    
    // Retrieve first list of CRT hits
    std::vector<crt::CRTHit> crtHitInList;
    art::Handle< std::vector<crt::CRTHit>> crtHitsInHandle;
    event.getByLabel(fCrtHitsIn_Label, crtHitsInHandle);
    //check to make sure the data we asked for is valid
    if(!crtHitsInHandle.isValid()){
      std::cout << "--> Run " << event.run() << ", subrun " << event.subRun()
		<< ", event " << event.event() << " has zero"
		<< " CRTHits " << " in module " << fCrtHitsIn_Label << std::endl;
      std::cout << std::endl;
      //add protection here
      event.put(std::move(CRTHitOutCol));
      return;
    }
    
    
    std::vector<crt::CRTHit> const& crtHitInList1(*crtHitsInHandle);
    crtHitInList.insert(crtHitInList.end(), crtHitInList1.begin(), crtHitInList1.end());
    if(fVerbose) std::cout<<"Number of CRT hits read in= "<<crtHitInList.size()<< 
		   " after first collection" << std::endl;
    
    

    
    for (size_t i = 0; i < crtHitInList.size(); i++){
      
      crt::CRTHit thisCrtHit = crtHitInList[i];
      std::vector<uint8_t> tfeb_id = thisCrtHit.feb_id; 
      double time1 = thisCrtHit.ts0_s;
      double time2 = thisCrtHit.ts0_s_corr;
      double time3 = thisCrtHit.ts0_ns;
      double time4 = thisCrtHit.ts0_ns_corr;
      double time5 = thisCrtHit.ts1_ns;
      
      int plane = thisCrtHit.plane;
      double x  = thisCrtHit.x_pos;
      double ex = thisCrtHit.x_err;
      double y  = thisCrtHit.y_pos;
      double ey = thisCrtHit.y_err;
      double z  = thisCrtHit.z_pos;
      double ez = thisCrtHit.z_err;

      std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap=thisCrtHit.pesmap;
      float pestot = thisCrtHit.peshit;      

      
      // only change/remove sim hits
      std::vector<std::pair<int,float>> test = tpesmap.find(tfeb_id[0])->second; 
      if (test.size()==2)  {
	time3 += thisCrtHit.ts1_ns;
      }


      // Create a corrected CRT hit
      crt::CRTHit crtHit = FillCrtHit(tfeb_id, tpesmap, pestot, time1,  time2,  time3,  time4,  time5, 
				      plane, x, ex,y,ey,z,ez );
      
      CRTHitOutCol->push_back(crtHit);
      nHits++;
            
    } // loop over hits
    
    event.put(std::move(CRTHitOutCol));

    if(fVerbose) std::cout<<"Number of CRT hits produced = "<<nHits<<std::endl;
      
    
  } // produce()
    



    void  CRTOverlayHotFix::endJob()
    {
      
    }
    
    crt::CRTHit  CRTOverlayHotFix::FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap, float peshit,double time1, double time2, double time3, double time4, double time5, 
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
      

      DEFINE_ART_MODULE( CRTOverlayHotFix)
      
  }// namespace crt

