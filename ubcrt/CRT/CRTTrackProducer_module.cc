////////////////////////////////////////////////////////////////////////
// Class:       CRTSimpleTrackProducer
// Module Type: producer
// File:        CRTSimpleTrackProducer_module.cc
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
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
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
#include <numeric>

namespace bernfebdaq {
  class CRTTrackProducer;
}

class bernfebdaq::CRTTrackProducer : public art::EDProducer {
public:
  explicit CRTTrackProducer(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  CRTTrackProducer(CRTTrackProducer const &) = delete;
  CRTTrackProducer(CRTTrackProducer &&) = delete;
  CRTTrackProducer & operator = (CRTTrackProducer const &) = delete;
  CRTTrackProducer & operator = (CRTTrackProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  //  art::ServiceHandle<art::TFileService> tfs;

  std::string  data_label_hits_;
  std::string  data_label_tzeros_;
  int track_method_type_;
  int store_track_;
  int verbose_ = 0;


  //quality plots
  //double track_time_ns = -1e18;
  //double track_time_s = -1e18;
  //double time_diff = 1e24;
  //double length = -1e18;
  //double theta = -1e18;
  //double phi = -1e18;

  //TH2F* hplavspla;
  //TH1F* hTlength;
  //TH1F* hTtime;
  //TH2F* hTlengthvsTime;
  //TH2F* hTlengthvsTimeAbs;
  //TProfile* hTlengthvsTimeAbs_prof;
  //TH1F* htheta;
  //TH1F* hphi;
  //TH1F* hts0_ns;
  //TH2F* hTvsH;

  //TH2F* HitDistBot;
  //TH2F* HitDistFT;
  //TH2F* HitDistPipe;
  //TH2F* HitDistTop;
 //quality plots                                                                                                     
 

};

void vmanip(std::vector<float> v, float* ave, float* rms);

struct CRTavehit{

  uint32_t ts0_ns;
  uint16_t ts0_ns_err;
  int32_t ts1_ns; 
  uint16_t ts1_ns_err;                                                        
  
  float x_pos;
  float x_err;
  float y_pos;
  float y_err;
  float z_pos;
  float z_err;
  float pe;

  int plane;
} tempah;


CRTavehit fillme(uint32_t i,uint16_t j,int32_t k,uint16_t l,float a,float b, float c,float d, float e, float f, float g, int p);

CRTavehit copyme(crt::CRTHit myhit);


crt::CRTTrack shcut(CRTavehit ppA,CRTavehit ppb, uint32_t time0s,uint16_t terr);



bernfebdaq::CRTTrackProducer::CRTTrackProducer(fhicl::ParameterSet const & p) : EDProducer{p}
{  
  // Initialize member data here.
  
  data_label_hits_ = p.get<std::string>("data_label_hits");
  data_label_tzeros_ = p.get<std::string>("data_label_tzeros");
  store_track_ =p.get<int>("store_track");
  verbose_ = p.get<int>("verbose");
  track_method_type_ = p.get<int>("track_method_type");
  
  
  // Call appropriate produces<>() functions here.
  if(store_track_ == 1) 
    produces< std::vector<crt::CRTTrack>   >();
}

void bernfebdaq::CRTTrackProducer::produce(art::Event & evt)
{

  //CRTTrack collection on this event                                                                         
  std::unique_ptr<std::vector<crt::CRTTrack> > CRTTrackCol(new std::vector<crt::CRTTrack>);

  art::Handle< std::vector<crt::CRTHit> > rawHandle;
  evt.getByLabel(data_label_hits_, rawHandle); //what is the product instance name? no BernZMQ
  
  //check to make sure the data we asked for is valid                                                                                                      
  if(!rawHandle.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTHits " << " in module " << data_label_hits_ << std::endl;
    std::cout << std::endl;

    if(store_track_ == 1)  evt.put(std::move(CRTTrackCol));

    return;
  }
  
  //get better access to the data               
  //  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle);


  art::Handle< std::vector<crt::CRTTzero> > rawHandletzero;
  evt.getByLabel(data_label_tzeros_, rawHandletzero); //what is the product instance name? no BernZMQ
  
  //check to make sure the data we asked for is valid                                                                                                      
  if(!rawHandletzero.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTTzeros " << " in module " << data_label_tzeros_ << std::endl;
    std::cout << std::endl;

    if(store_track_ == 1) evt.put(std::move(CRTTrackCol));

    return;
  }
  
  //get better access to the data               
  //  std::vector<crt::CRTTzero> const& CRTTCollection(*rawHandletzero);

  //CRTTzero collection on this event                                                                         
  std::unique_ptr<std::vector<crt::CRTTzero> > CRTTzeroCol(new std::vector<crt::CRTTzero>);
  
    std::vector<art::Ptr<crt::CRTTzero> > tzerolist;

  if (evt.getByLabel(data_label_tzeros_,rawHandletzero))
    art::fill_ptr_vector(tzerolist, rawHandletzero);

  art::FindManyP<crt::CRTHit> fmht(rawHandletzero, evt, data_label_tzeros_);


  //loop over tzeros
  for(size_t tzIter = 0; tzIter < tzerolist.size(); ++tzIter){   
    
    //count planes with hits for this tzero
    int np =0 ; int ipflag[4] = {};
    int tothits =0;
    for (int ip=0;ip<4;++ip) {
      if (tzerolist[tzIter]->nhits[ip]>0)  { np++; ipflag[ip]=1; tothits+=tzerolist[tzIter]->nhits[ip];}  
    }

    if (np>1) {
    std::vector<art::Ptr<crt::CRTHit> > hitlist=fmht.at(tzIter);
      if (track_method_type_==1) {

	double time_s_A = hitlist[0]->ts0_s;
	// find pairs of hits in different planes
	for (size_t ah = 0; ah< hitlist.size()-1; ++ah){	
	  crt::CRTHit temphit=*hitlist[ah];
	  CRTavehit Ahit = copyme(temphit);
	  int planeA = hitlist[ah]->plane%10;
	  for (size_t bh = ah+1; bh< hitlist.size(); ++bh){	
	    int planeB = hitlist[bh]->plane%10;
	    if (planeB!=planeA) {  // make a track	       
	      temphit=*hitlist[bh];
	      CRTavehit Bhit = copyme(temphit);
	      crt::CRTTrack CRTcanTrack=shcut(Ahit,Bhit,time_s_A,0);
	      CRTTrackCol->emplace_back(CRTcanTrack);
	    }
	  }
	}
      }
      else if ((track_method_type_==2) || (track_method_type_==3 && np==2 && tothits==2)) {	
	//loop over hits and get average x,y,z,pe for each plane
	std::vector<float> thittime0[4];
	std::vector<float> thittime1[4];
	std::vector<float> tx[4];
	std::vector<float> ty[4];
	std::vector<float> tz[4];
	std::vector<float> pe[4];
	
	double time_s_A = hitlist[0]->ts0_s;
	//      double time_s_err = hitlist[0]->ts0_s_err;
	double time_s_err = 0.;
	double time1_ns_A = hitlist[0]->ts1_ns;
	double time0_ns_A = hitlist[0]->ts0_ns;
	
	//loop over hits for this tzero, sort by plane
	for (size_t ah = 0; ah< hitlist.size(); ++ah){	
	  int ip = hitlist[ah]->plane%10;       
	  thittime0[ip].push_back(hitlist[ah]->ts0_ns-time0_ns_A);
	  thittime1[ip].push_back(hitlist[ah]->ts1_ns-time1_ns_A);
	  tx[ip].push_back(hitlist[ah]->x_pos);
	  ty[ip].push_back(hitlist[ah]->y_pos);
	  tz[ip].push_back(hitlist[ah]->z_pos);
	  pe[ip].push_back(hitlist[ah]->peshit);	
	} // loop over hits
	
	  // now take averages if there are multiple hits in the same plane at the same time
	uint planeA,planeB,planeC,planeD;
	float totpe=0.0;
	float avet1=0.0; float rmst1 =0.0; 
	float avet0=0.0; float rmst0 =0.0; 
	float avex=0.0; float rmsx =0.0; 
	float avey=0.0; float rmsy =0.0; 
	float avez=0.0; float rmsz =0.0; 
	uint32_t at0;        int32_t at1;   uint16_t rt0,rt1;
	planeA=-1;      planeB=-1;      planeC=-1;      planeD=-1;
	int ip=0;
	while (ip<4 && ipflag[ip]==0) ip++;
	planeA=ip;
	ip++;
	while (ip<4 && ipflag[ip]==0) ip++;
	planeB=ip;
	
	//First track A-B
	///average over hits in plane A
	vmanip(thittime0[planeA],&avet0,&rmst0);
	vmanip(thittime1[planeA],&avet1,&rmst1);
	at0 = (uint32_t)(avet0+time0_ns_A); rt0 = (uint16_t)rmst0;   
	at1 = (int32_t)(avet1+time1_ns_A); rt1 = (uint16_t)rmst1;
	vmanip(tx[planeA],&avex,&rmsx);
	vmanip(ty[planeA],&avey,&rmsy);
	vmanip(tz[planeA],&avez,&rmsz);
	totpe=std::accumulate(pe[planeA].begin(), pe[planeA].end(), 0.0);
	CRTavehit pA = fillme(at0,rt0,at1,rt1,avex,rmsx,avey,rmsy,avez,rmsz,totpe,planeA);      
	//  average over hits in plane B
	vmanip(thittime0[planeB],&avet0,&rmst0);
	vmanip(thittime1[planeB],&avet1,&rmst1);
	at0 = (uint32_t)(avet0+time0_ns_A); rt0 = (uint16_t)rmst0;   
	at1 = (int32_t)(avet1+time1_ns_A); rt1 = (uint16_t)rmst1;
	vmanip(tx[planeB],&avex,&rmsx);
	vmanip(ty[planeB],&avey,&rmsy);
	vmanip(tz[planeB],&avez,&rmsz);
	totpe=std::accumulate(pe[planeB].begin(), pe[planeB].end(), 0.0);
	CRTavehit pB = fillme(at0,rt0,at1,rt1,avex,rmsx,avey,rmsy,avez,rmsz,totpe,planeB);      
	crt::CRTTrack CRTcanTrack=shcut(pA,pB,time_s_A,time_s_err);
	CRTTrackCol->emplace_back(CRTcanTrack);
	//second and third track
	if (np>2 && ip<3) {
	  ip++;
	  while (ip<4 && ipflag[ip]==0) ip++;
	  planeC=ip;
	  //	std::cout << "plane C is " << planeC << std::endl;
	  vmanip(thittime0[planeC],&avet0,&rmst0);
	  vmanip(thittime1[planeC],&avet1,&rmst1);
	  at0 = (uint32_t)(avet0+time0_ns_A); rt0 = (uint16_t)rmst0;   
	  at1 = (int32_t)(avet1+time1_ns_A); rt1 = (uint16_t)rmst1;
	  vmanip(tx[planeC],&avex,&rmsx);
	  vmanip(ty[planeC],&avey,&rmsy);
	  vmanip(tz[planeC],&avez,&rmsz);
	  totpe=std::accumulate(pe[planeC].begin(), pe[planeC].end(), 0.0);
	  CRTavehit pC = fillme(at0,rt0,at1,rt1,avex,rmsx,avey,rmsy,avez,rmsz,totpe,planeC);
	  CRTcanTrack=shcut(pA,pC,time_s_A,time_s_err);
	  CRTTrackCol->emplace_back(CRTcanTrack);
	  CRTcanTrack=shcut(pB,pC,time_s_A,time_s_err);
	  CRTTrackCol->emplace_back(CRTcanTrack);
	  // last 3 tracks
	  if (np==4) {
	    planeD=3;
	    //	  std::cout << "plane D is " << planeD << std::endl;
	    vmanip(thittime0[planeD],&avet0,&rmst0);
	    vmanip(thittime1[planeD],&avet1,&rmst1);
	    at0 = (uint32_t)(avet0+time0_ns_A); rt0 = (uint16_t)rmst0;   
	    at1 = (int32_t)(avet1+time1_ns_A); rt1 = (uint16_t)rmst1;
	    vmanip(tx[planeD],&avex,&rmsx);
	    vmanip(ty[planeD],&avey,&rmsy);
	    vmanip(tz[planeD],&avez,&rmsz);
	    totpe=std::accumulate(pe[planeD].begin(), pe[planeD].end(), 0.0);
	    CRTavehit pD =fillme(at0,rt0,at1,rt1,avex,rmsx,avey,rmsy,avez,rmsz,totpe,planeD);	  
	    CRTcanTrack=shcut(pA,pD,time_s_A,time_s_err);
	    CRTTrackCol->emplace_back(CRTcanTrack);
	    CRTcanTrack=shcut(pB,pD,time_s_A,time_s_err);
	    CRTTrackCol->emplace_back(CRTcanTrack);
	    CRTcanTrack=shcut(pC,pD,time_s_A,time_s_err);
	    CRTTrackCol->emplace_back(CRTcanTrack);
	  }
	}
      }
	
    }// if at least two planes with hits
    
  }// loop over tzeros
  
  //store track collection into event
  if(store_track_ == 1)
    evt.put(std::move(CRTTrackCol));
  
}

void bernfebdaq::CRTTrackProducer::beginJob()
{
  

}

void bernfebdaq::CRTTrackProducer::endJob()
{
  // Implementation of optional member function here.
}



void vmanip(std::vector<float> v, float* ave, float* rms)
{
  *ave=0.0; *rms =0.0;
  if (v.size()>0) {

    /*
    int np=v.size();
    for (int i=0;i<np;++i) {
      std::cout << v[i] << " " ;
    }
    std::cout << std::endl;
    */

    //  find the mean and *rms of all the vector elements
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    *ave=mean;
    
    if (v.size()>1) {
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
    *rms=stdev;
    }

  }

  //    std::cout << "inside vmanip " << *ave << " " << *rms << std::endl;
}


CRTavehit fillme(uint32_t i,uint16_t j,int32_t k,uint16_t l,float a,float b, 
		 float c, float d, float e,float f,float g,int p)
{

  CRTavehit h;
  h.ts0_ns=i;
  h.ts0_ns_err=j;
  h.ts1_ns=k; 
  h.ts1_ns_err=l;                                                        
  
  h.x_pos=a;
  h.x_err=b;
  h.y_pos=c;
  h.y_err=d;
  h.z_pos=e;
  h.z_err=f;
  h.pe=g;
  h.plane=p;
  return(h);


} 


CRTavehit copyme(crt::CRTHit myhit)
{

  CRTavehit h;
  h.ts0_ns=myhit.ts0_ns;
  h.ts0_ns_err=0;
  h.ts1_ns=myhit.ts1_ns;; 
  h.ts1_ns_err=0;       
  h.x_pos=myhit.x_pos;
  h.x_err=myhit.x_err;
  h.y_pos=myhit.y_pos;
  h.y_err=myhit.y_err;
  h.z_pos=myhit.z_pos;
  h.z_err=myhit.z_err;
  h.pe=myhit.peshit;
  h.plane=myhit.plane%10;
  return(h);


} 



crt::CRTTrack shcut(CRTavehit ppA,CRTavehit ppB,uint32_t time0s,uint16_t terr)
{

  crt::CRTTrack newtr;
  newtr.ts0_s=time0s;
  newtr.ts0_s_err=terr;
  newtr.ts0_ns_h1=ppA.ts0_ns;
  newtr.ts0_ns_err_h1=ppA.ts0_ns_err;
  newtr.ts0_ns_h2=ppB.ts0_ns;
  newtr.ts0_ns_err_h2=ppB.ts0_ns_err;
  newtr.ts0_ns=(uint32_t)(0.5*(ppA.ts0_ns+ppB.ts0_ns));
  newtr.ts0_ns_err=(uint16_t)(0.5*sqrt(ppA.ts0_ns_err*ppA.ts0_ns_err+ppB.ts0_ns_err*ppB.ts0_ns_err));
  newtr.ts1_ns=(int32_t)(0.5*(ppA.ts1_ns+ppB.ts1_ns));
  newtr.ts1_ns_err=(uint16_t)(0.5*sqrt(ppA.ts0_ns_err*ppA.ts0_ns_err+ppB.ts0_ns_err*ppB.ts0_ns_err));
  newtr.peshit=ppA.pe+ppB.pe;
  newtr.x1_pos=ppA.x_pos;
  newtr.x1_err=ppA.x_err;
  newtr.y1_pos=ppA.y_pos;
  newtr.y1_err=ppA.y_err;
  newtr.z1_pos=ppA.z_pos;
  newtr.z1_err=ppA.z_err;
  newtr.x2_pos=ppB.x_pos;
  newtr.x2_err=ppB.x_err;
  newtr.y2_pos=ppB.y_pos;
  newtr.y2_err=ppB.y_err;
  newtr.z2_pos=ppB.z_pos;
  newtr.z2_err=ppB.z_err;
  float deltax =  ppA.x_pos-ppB.x_pos;
  float deltay =  ppA.y_pos-ppB.y_pos;
  float deltaz =  ppA.z_pos-ppB.z_pos;
  newtr.length=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
  newtr.thetaxy=atan2(deltax,deltay);
  newtr.phizy=atan2(deltaz,deltay);
  newtr.plane1=ppA.plane%10;
  newtr.plane2=ppB.plane%10;
  return(newtr);

}


DEFINE_ART_MODULE(bernfebdaq::CRTTrackProducer)
