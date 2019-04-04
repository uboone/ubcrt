////////////////////////////////////////////////////////////////////////
// Class:       CRTTzeroProducer
// Module Type: producer
// File:        CRTTzeroProducer_module.cc
// Description: Module for constructiong over-simplified CRT tracks.
// Copied from CRTTrackProducer by David Lorca Galindo 
//  Edited by Michelle Stancari April 3, 2018
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CRTBernFEBDAQCore/Overlays/BernZMQFragment.hh"
#include <artdaq-core/Data/Fragment.hh>

#include "art_root_io/TFileService.h"

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTTrack.hh"
#include "ubobj/CRT/CRTTzero.hh"
#include "ubcrt/CRT/CRTAuxFunctions.hh"

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
#include <cmath> 
#include <memory>

// namespace bernfebdaq {
//   class CRTTzeroProducer;
// }

class CRTTzeroProducer : public art::EDProducer {
public:
  explicit CRTTzeroProducer(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  CRTTzeroProducer(CRTTzeroProducer const &) = delete;
  CRTTzeroProducer(CRTTzeroProducer &&) = delete;
  CRTTzeroProducer & operator = (CRTTzeroProducer const &) = delete;
  CRTTzeroProducer & operator = (CRTTzeroProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  //  art::ServiceHandle<art::TFileService> tfs;

  std::string  data_label_;
  double max_time_difference_ ;//max time for coincidence 
  int store_tzero_;
  int verbose_ = 0;

};

void vmanip(std::vector<double> v, double* ave, double* rms);
void set_def(crt::CRTTzero tz);


CRTTzeroProducer::CRTTzeroProducer(fhicl::ParameterSet const & p)
  :
  // Initialize member data here.
  data_label_(p.get<std::string>("data_label")),
  max_time_difference_(p.get<double>("max_time_difference")),
  store_tzero_(p.get<int>("store_tzero")),
  verbose_(p.get<int>("verbose"))
  
{
  // Call appropriate produces<>() functions here.
  if(store_tzero_ == 1)  {
    produces< std::vector<crt::CRTTzero>   >();
    produces<art::Assns<crt::CRTTzero, crt::CRTHit> >();
  }
}

void CRTTzeroProducer::produce(art::Event & evt)
{
  // Implementation of required member function here.
  
  art::Handle< std::vector<crt::CRTHit> > rawHandle;
  evt.getByLabel(data_label_, rawHandle);   
  //check to make sure the data we asked for is valid                                                                                                      
  if(!rawHandle.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTHits " << " in module " << data_label_ << std::endl;
    std::cout << std::endl;
    return;
  }
  
  //get better access to the data               
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle);

  //CRTTzero collection on this event                                                              
  std::unique_ptr<std::vector<crt::CRTTzero> > CRTTzeroCol(new std::vector<crt::CRTTzero>);
  

  // Output collections  
  std::unique_ptr<art::Assns<crt::CRTTzero, crt::CRTHit>> outputHits(new art::Assns<crt::CRTTzero, crt::CRTHit>);
  //  auto outputHits    = std::make_unique<art::Assns<crt::CRTTzero, crt::CRTHit>>();
  // need later version of art (later than v2_05_01) to use PtrMaker.
  // do it the old-fashioned way instead
  //      art::PtrMaker<crt::CRTHit> hitPtrMaker(evt, rawhandle.id());
  //  art::PtrMaker<crt::CRTTzero> tzeroPtrMaker(evt);
  art::PtrMaker<crt::CRTHit> hitPtrMaker(evt, rawHandle.id());
  art::PtrMaker<crt::CRTTzero> tzeroPtrMaker(evt);

 

  int N_CRTHits = CRTHitCollection.size();
  //  std::cout << "number of crt hits " << N_CRTHits << std::endl;
  int iflag[1000] = {};

  uint planeA, planeB;


  for(int  i = 0; i < N_CRTHits; i++) {//A 
        
    if (iflag[i]==0) {  // new tzero
      //temporary hit collection for each tzero
      std::vector<art::Ptr<crt::CRTHit>> CRTHitCol;
      crt::CRTHit CRTHiteventA = CRTHitCollection[i];
      art::Ptr<crt::CRTHit> hptr = hitPtrMaker(i);
      CRTHitCol.push_back(hptr);
      double time_s_A = CRTHiteventA.ts0_s;
      double time_ns_A = CRTHiteventA.ts1_ns;
      iflag[i]=1;
      // create and initialize, ugly code :(
      crt::CRTTzero CRTcanTzero;
      CRTcanTzero.ts0_ns=0;
      CRTcanTzero.ts1_ns=0;
      for (int i=0;i<4;++i) {
	CRTcanTzero.nhits[i]=0;
	CRTcanTzero.pes[i]=0;
      }
      // 
      CRTcanTzero.ts0_s=CRTHiteventA.ts0_s;
      CRTcanTzero.ts0_s_err=0;
      planeA = CRTHiteventA.plane;
      CRTcanTzero.nhits[planeA]=1;      
      int icount=1;
      CRTcanTzero.pes[planeA]=CRTHiteventA.peshit;
      for(int j = i+1; j < N_CRTHits; j++) {//B
	if (iflag[j]==0) {
	  crt::CRTHit CRTHiteventB = CRTHitCollection[j];
	  //look for coincidences
	  double time_s_B = CRTHiteventB.ts0_s;
	  double time_ns_B = CRTHiteventB.ts1_ns;
	  double time_diff = time_ns_B - time_ns_A;
	  if( (time_s_A == time_s_B) && (abs(time_diff)<max_time_difference_)  ){//D
	    art::Ptr<crt::CRTHit> hptr = hitPtrMaker(j);
	    CRTHitCol.push_back(hptr);
	    planeB = CRTHiteventB.plane; 
	    CRTcanTzero.nhits[planeB]+=1;
	    CRTcanTzero.pes[planeB]+=CRTHiteventB.peshit;
	    iflag[j]=1;
	    CRTcanTzero.ts1_ns+=(int)(time_diff);
	    CRTcanTzero.ts0_ns+=(CRTHiteventB.ts0_ns-CRTHiteventA.ts0_ns);
	    icount++;
	  }
	}
      }      // done with this tzero


      // Make a tzero data product
      CRTcanTzero.ts1_ns/=icount; 
      CRTcanTzero.ts1_ns+=(int)time_ns_A;
      CRTcanTzero.ts0_ns/=icount;
      CRTcanTzero.ts0_ns+=(int)CRTHiteventA.ts0_ns;
      CRTcanTzero.ts1_ns_err=0.;
      CRTcanTzero.ts0_ns_err=0.;
      CRTTzeroCol->push_back(CRTcanTzero);
      //associate hits to this Tzero
      art::Ptr<crt::CRTTzero> aptz = tzeroPtrMaker(CRTTzeroCol->size()-1);
      util::CreateAssn(*this,evt,aptz,CRTHitCol,*outputHits);

    }//B
  }//A
  
  
  //store tzero collection into event
  if(store_tzero_ == 1) {
    evt.put(std::move(CRTTzeroCol));
    evt.put(std::move(outputHits));
  }
  
}

void CRTTzeroProducer::beginJob()
{
  

}

void CRTTzeroProducer::endJob()
{
  // Implementation of optional member function here.
}




DEFINE_ART_MODULE(CRTTzeroProducer)


