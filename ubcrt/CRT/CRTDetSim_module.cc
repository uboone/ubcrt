#include "ubobj/CRT/CRTSimData.hh"
#include "ubcrt/CRT/CRTDetSim.hh"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()
#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

#include <vector>
#include <memory>
#include <sstream>
#include <string>
#include <math.h>

namespace crt{
  /*
  const double striplength[73] = {
    346.0, 346.0, 346.0, 259.6, 259.6, 259.6, 259.6, 227.0, 227.0, 346.0, 
    346.0, 346.0, 346.0, 346.0, 346.0, 346.0, 403.8, 403.8, 403.8, 403.8, 
    403.8, 403.8, 259.6, 259.6, 259.6, 259.6, 259.6, 259.6, 259.6, 259.6, 
    259.6, 259.6, 259.6, 259.6, 259.6, 259.6, 396.2, 396.2, 396.2, 227.0, 
    227.0, 227.0, 227.0, 227.0, 227.0, 227.0, 227.0, 227.0, 227.0, 360.0, 
    360.0, 360.0, 360.0, 360.0, 360.0, 360.0, 180.0, 180.0, 180.0, 180.0, 
    180.0, 180.0, 365.0, 365.0, 365.0, 365.0, 365.0, 365.0, 180.0, 180.0, 
    365.0, 365.0, 365.0};*/
    const short mod2plane[73] = {
      0,0,0,0,0,0,0,0,0,1,  //0-9
      1,1,1,1,1,1,1,1,1,1,  //10-19
      1,1,2,2,2,2,2,2,2,2,  //20-29
      2,2,2,2,2,2,2,2,2,2,  //30-39
      2,2,2,2,2,2,2,2,2,3,  //40-49
      3,3,3,3,3,3,3,3,3,3,  //50-59  
      3,3,3,3,3,3,3,3,3,3,  //60-69  
      3,3,3};               //70-72 
  /* 
    const short mod2feb[73] = {
      24,23,22,17,14,18,19,12,11,52,  //0-9
      31,29,28,27,26,30,61,59,57,60,                  //10-19
      58,56,32,38,36,35,34,33,37,45,                  //20-29
      44,43,42,41,40,39,55,54,53,51,                  //30-39
      49,47,21,16,50,48,46,20,15,107,                 //40-49
      106,105,109,108,112,111,115,123,124,126,        //50-59
      125,129,115,114,113,116,119,120,127,128,        //60-69
      117,121,118};                                   //70-72
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


  CRTDetSim::CRTDetSim(const fhicl::ParameterSet& pSet)
  {
    art::ServiceHandle<rndm::NuRandomService> Seeds;
    Seeds->createEngine(*this, "HepJamesRandom", "crt", pSet, "Seed");
    this->reconfigure(pSet);
    produces< std::vector<CRTSimData> >();

  }

  CRTDetSim::~CRTDetSim()
  {

  }

  void CRTDetSim::reconfigure(fhicl::ParameterSet const & pSet) {
    fG4ModuleLabel = pSet.get<std::string>("G4ModuleLabel","g4");
    fGlobalT0Offset = pSet.get<double>("GlobalT0Offset",0.0);
    fTDelayNorm = pSet.get<double>("TDelayNorm",4126.);
    fTDelayShift = pSet.get<double>("TDelayShift",-300.);
    fTDelaySigma = pSet.get<double>("TDelaySigma",90.);
    fTDelayOffset = pSet.get<double>("TDelayOffset",-1.5);
    fTDelayRMSGausNorm = pSet.get<double>("TDelayRMSGausNorm",2.09);
    fTDelayRMSGausShift = pSet.get<double>("TDelayRMSGausShift",7.24);
    fTDelayRMSGausSigma = pSet.get<double>("TDelayRMSGausSigma",170.);
    fTDelayRMSExpNorm = pSet.get<double>("TDelayRMSExpNorm",1.65);
    fTDelayRMSExpShift = pSet.get<double>("TDelayRMSExpShift",75.6);
    fTDelayRMSExpScale = pSet.get<double>("TDelayRMSExpScale",79.35);
    fPropDelay = pSet.get<double>("PropDelay",0.0061);
    fPropDelayError = pSet.get<double>("PropDelayError",0.007);
    fTResInterpolator = pSet.get<double>("TResInterpolator",1.268);
    fNpeScaleNorm = pSet.get<double>("NpeScaleNorm",5261000.);
    fNpeScaleShift = pSet.get<double>("NpeScaleShift",-1085.);
    fQ0 = pSet.get<double>("Q0",0.00175);
    fQPed = pSet.get<double>("QPed",63.6);
    fQSlope = pSet.get<double>("QSlope",131.9);
    fQRMS = pSet.get<double>("QRMS",15.);
    fQThreshold = pSet.get<double>("QThreshold",100.0);
    fAbsLenEff = pSet.get<double>("AbsLenEff",8.5);
    fStripCoincidenceWindow = pSet.get<double>("StripCoincidenceWindow", 30.);
    fTaggerPlaneCoincidenceWindow = pSet.get<double>("TaggerPlaneCoincidenceWindow", 5.);
    fSipmTimeResponse = pSet.get<double>("SipmTimeResponse",2.);
    fUseEdep = pSet.get<bool>("UseEdep",true);
    fverbose = pSet.get<bool>("verbose",0);

    if (fQThreshold<fQPed) { 
      fQThreshold=fQPed;  
      std::cout << " Threshold cannot be lower than pedestal value.  Setting threshold at pedestal" << std::endl;
    }
  }

  double CRTDetSim::getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                           detinfo::ElecClock& clock,
                                           float t0, float npeMean, float r) {
    // Hit timing, with smearing and NPE dependence
    double tDelayMean = \
      fTDelayNorm *
        exp(-0.5 * pow((npeMean - fTDelayShift) / fTDelaySigma, 2)) +
      fTDelayOffset;

    double tDelayRMS = \
      fTDelayRMSGausNorm *
        exp(-pow(npeMean - fTDelayRMSGausShift, 2) / fTDelayRMSGausSigma) +
      fTDelayRMSExpNorm *
        exp(-(npeMean - fTDelayRMSExpShift) / fTDelayRMSExpScale);

    double tDelay = CLHEP::RandGauss::shoot(engine, tDelayMean, tDelayRMS);

    // Time resolution of the interpolator
    tDelay += CLHEP::RandGauss::shoot(engine, 0, fTResInterpolator);

    // Propagation time
    double tProp = CLHEP::RandGauss::shoot(fPropDelay, fPropDelayError) * r;
    double t = t0 + tProp + tDelay;

    // Get clock ticks
    clock.SetTime(t/1000.);  // SetTime takes microseconds

    mf::LogInfo("CRT")
    << "CRT TIMING: t0=" << t0
    << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS
    << ", tDelay=" << tDelay << ", tDelay(interp)="
    << tDelay << ", tProp=" << tProp << ", t=" << t << ", ticks=" << clock.Ticks() << "\n";


    return clock.Ticks();
  }


  struct Tagger {
    std::vector<std::pair<unsigned, uint32_t> > planesHit;
    std::vector<crt::CRTSimData> data;
  };

  void CRTDetSim::produce(art::Event& evt)
  {
   // A list of hit taggers, before any coincidence requirement
    std::map<std::string, Tagger> taggers;

    std::unique_ptr<std::vector<crt::CRTSimData> > crtHits(
        new std::vector<crt::CRTSimData>);

    art::ServiceHandle<detinfo::DetectorClocksService> detClocks;
    detinfo::ElecClock trigClock = detClocks->provider()->TriggerClock();

    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine("crt");

    // Handle for (truth) AuxDetSimChannels
    art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
    evt.getByLabel(fG4ModuleLabel, channels);

    //access geometry
    art::ServiceHandle<geo::AuxDetGeometry> geoService;
    const geo::AuxDetGeometry* geometry = &*geoService;
    const geo::AuxDetGeometryCore* geoServiceProvider = geometry->GetProviderPtr();

    // Loop through truth AD channels
    for (auto& adsc : *channels) {
      const geo::AuxDetGeo& adGeo = geoServiceProvider->AuxDet(adsc.AuxDetID());
      const geo::AuxDetSensitiveGeo& adsGeo = adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());
 
    // Return the vector of IDEs
      std::vector<sim::AuxDetIDE> ides = adsc.AuxDetIDEs();
      std::sort(ides.begin(), ides.end(), 
		[](const sim::AuxDetIDE & a, const sim::AuxDetIDE & b) -> bool{
		  return ((a.entryT + a.exitT)/2) < ((b.entryT + b.exitT)/2);
		});

      // Simulate the CRT response for each hit
      for (size_t ide_i = 0; ide_i < ides.size(); ide_i++) {
	sim::AuxDetIDE ide = ides[ide_i];
   	int trackID = ide.trackID;
           // Get the hit position in strip's local coordinates
        double x = (ide.entryX + ide.exitX) / 2;
        double y = (ide.entryY + ide.exitY) / 2;
        double z = (ide.entryZ + ide.exitZ) / 2;

        double tTrue = (ide.entryT + ide.exitT) / 2 + fGlobalT0Offset;
        double eDep = ide.energyDeposited;	
	double maxEdep = eDep;	


	//ADD UP HITS AT THE SAME TIME - 2NS DIFF IS A GUESS -VERY APPROXIMATE
	if(ide_i< ides.size() - 1){
	  while(ide_i < ides.size() - 1 && std::abs(tTrue-((ides[ide_i+1].entryT + ides[ide_i+1].exitT) / 2 + fGlobalT0Offset)) < fSipmTimeResponse){
	    ide_i++;
	    x = (x + ((ides[ide_i].entryX + ides[ide_i].exitX) / 2))/2;
	    y = (y + ((ides[ide_i].entryY + ides[ide_i].exitY) / 2))/2;
	    z = (z + ((ides[ide_i].entryZ + ides[ide_i].exitZ) / 2))/2;
	    eDep += ides[ide_i].energyDeposited;
	    tTrue = (tTrue + ((ides[ide_i].entryT + ides[ide_i].exitT) / 2))/2;
	    if(ides[ide_i].energyDeposited > maxEdep && ides[ide_i].trackID > 0){
	      trackID = ides[ide_i].trackID;
	      maxEdep = ides[ide_i].energyDeposited;
	    }
	    if(trackID < 0 && (ides[ide_i].energyDeposited > maxEdep || ides[ide_i].trackID > 0)){
	      trackID = ides[ide_i].trackID;
            maxEdep = ides[ide_i].energyDeposited;
	    }
	  }
	}

	std::string name = adGeo.TotalVolume()->GetName();
	int module,strip;
	sscanf(name.c_str(),"volAuxDet_Module_%d_strip_%d",&module,&strip);
	if (fverbose) std::cout << " channel " << adsc.AuxDetID() << " strip " << strip << 
			" module " << module << " name " << name << std::endl;
	double world[3] = {x, y, z};
	double svHitPosLocal[3];
	adsGeo.WorldToLocal(world, svHitPosLocal);	

        // Distance to the readout end 
        // FIXME: FOR NOW ASSUME ALL THE SAME DIRECTION (SiPM always at smaller coordinate value)
	// FIXME: the y coordinate is not always along the strip length in uB like it is for SBND.
	double distToReadout= abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);
	if (mod2orient[module]==0)   distToReadout = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);
	else if (mod2orient[module]==2)   distToReadout = abs(-adsGeo.HalfLength() - svHitPosLocal[2]);
      
	// The expected number of PE, using a quadratic model for the distance
	// dependence, and scaling linearly with deposited energy.
	//  UseEdep flag let's us do a binary response
	double qr = fUseEdep ? 1.0 * eDep / fQ0 : 1.0;
	
	double npeExpected = (fNpeScaleNorm / pow(distToReadout - fNpeScaleShift, 2) * qr);
	
	// Put PE on channels weighted by distance
	double d0 = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);  // L
	double d1 = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]);  // R
	if (mod2orient[module]==0) {
	  if (mod2plane[module]==0  || mod2plane[module]==3 ) {
	    d0 = abs(-adsGeo.HalfLength() - svHitPosLocal[2]);  // L
	    d1 = abs( adsGeo.HalfLength() - svHitPosLocal[2]);  // R
	  }
	  else {
	    d0 = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);  // L
	    d1 = abs( adsGeo.HalfHeight() - svHitPosLocal[1]);  // R	    
	  }
	}
	else if (mod2orient[module]==z) {
	  if (mod2plane[module]==0  || mod2plane[module]==3 ) {
	    d0 = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);  // L
	    d1 = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]);  // R
	  }
	  else {
	    d0 = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);  // L
	    d1 = abs( adsGeo.HalfHeight() - svHitPosLocal[1]);  // R	    
	  }
	}
	double abs0 = exp(-d0 / fAbsLenEff);
	double abs1 = exp(-d1 / fAbsLenEff);
	double npeExp0 = npeExpected * abs0 / (abs0 + abs1);
	double npeExp1 = npeExpected * abs1 / (abs0 + abs1);
	
        // Observed PE
        long npe0 = CLHEP::RandPoisson::shoot(engine, npeExp0);
        long npe1 = CLHEP::RandPoisson::shoot(engine, npeExp1);

        // Time relative to trigger
        uint32_t t0 = getChannelTriggerTicks(engine, trigClock, tTrue, npe0, distToReadout);
        uint32_t t1 = getChannelTriggerTicks(engine, trigClock, tTrue, npe1, distToReadout);

        // Time relative to PPS: Random for now
        uint32_t ppsTicks = CLHEP::RandFlat::shootInt(engine, trigClock.Frequency() * 1e6);

        // SiPM and ADC response: Npe to ADC counts
        short q0 = CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe0, fQRMS * sqrt(npe0));
        short q1 = CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe1, fQRMS * sqrt(npe1));

        uint32_t channel0ID = adsc.AuxDetID()*2;
        uint32_t channel1ID = adsc.AuxDetID()*2+1;

	if (fverbose) 
	  std::cout << adsc.AuxDetID() << " " << ide_i << " " << channel0ID << " " << q0 << 
	    " " << t0 <<
	    " " << channel1ID << " " << q1 << " " << t1 <<  std::endl;

      // Apply ADC threshold and strip-level coincidence (both fibers fire)
	if (q0 > fQThreshold &&
	    q1 > fQThreshold &&
	    util::absDiff(t0, t1) < fStripCoincidenceWindow) {
	   crtHits->push_back(CRTSimData(channel0ID, t0, ppsTicks, q0, trackID));
	   crtHits->push_back(CRTSimData(channel1ID, t1, ppsTicks, q1, trackID));
	   if (fverbose) std::cout << "both written to event q0: " << q0 << " q1 " << q1 <<std::endl;
	}
      } //loop over ides
    }// loop over channels

    // trigger commented out because it is duplicated in the hit reconstruction code CRTSimHit_module.cc
    // this code does not work because it is SBND specific, but could be modified easily
    //
    // Apply coincidence trigger requirement
    // Logic: require at least one hit in each perpendicular plane. 
    
    /*
      std::unique_ptr<std::vector<sbnd::crt::CRTData> > triggeredCRTHits(
      new std::vector<sbnd::crt::CRTData>);
      // Logic: For normal taggers, require at least one hit in each perpendicular
      // plane. For the bottom tagger, any hit triggers read out.
      for (auto trg : taggers) {
      bool trigger = false;
      // Loop over pairs of hits
      for (auto t1 : trg.second.planesHit) {
      for (auto t2 : trg.second.planesHit) {
      if (t1 == t2) {
      continue;
      }
      // Two hits on different planes with proximal t0 times
      if (t1.first != t2.first && std::abs(t1.second - t2.second) < fTaggerPlaneCoincidenceWindow) {
      trigger = true;
      break;
      }
      }
    }
    if (trigger || trg.first.find("TaggerBot") != std::string::npos) {
    // Write out all hits on a tagger when there is any coincidence
    for (auto d : trg.second.data) {
    triggeredCRTHits->push_back(d);
    }
    }
    } // end loop over hits
    mf::LogInfo("CRT") << "CRT TRIGGERED HITS: " << triggeredCRTHits->size() << "\n";
    e.put(std::move(triggeredCRTHits));
    */
  

    evt.put(std::move(crtHits));
  }
  

  DEFINE_ART_MODULE(CRTDetSim)
}  //namespace CRT
