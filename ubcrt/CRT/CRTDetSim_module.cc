/**
 * \class CRTDetSim
 *
 * \ingroup crt
 *
 * \brief Provides CRTData from simulations
 *
 * Converts IDEs from largeant (or whichever producer) to
 * CRTData. This is meant to mimic the physical detector as much as
 * possible.
 *
 *
 * \author $Author: Kevin Wierman<kevin.wierman@pnnl.gov>
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2016/12/12 $
 *
 * Contact: kevin.wierman@pnnl.gov
 *
 * Created on: Tuesday, December 13, 2016
 *
**/

#include "ubobj/CRT/CRTSimData.hh"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

#include <vector>
#include <memory>
#include <sstream>
#include <string>
#include <math.h>

namespace crt{
  class CRTDetSim :  public art:: EDProducer{
  public:
    explicit CRTDetSim(const fhicl::ParameterSet&);

  private:
    void produce (art::Event&) override;

    /// Name of the producer of the IDEs
    std::string fG4ModuleLabel;
    double fGlobalT0Offset;  //!< Time delay fit: Gaussian normalization
    double fTDelayNorm;  //!< Time delay fit: Gaussian normalization
    double fTDelayShift;  //!< Time delay fit: Gaussian x shift
    double fTDelaySigma;  //!< Time delay fit: Gaussian width
    double fTDelayOffset;  //!< Time delay fit: Gaussian baseline offset
    double fTDelayRMSGausNorm;  //!< Time delay RMS fit: Gaussian normalization
    double fTDelayRMSGausShift;  //!< Time delay RMS fit: Gaussian x shift
    double fTDelayRMSGausSigma;  //!< Time delay RMS fit: Gaussian width
    double fTDelayRMSExpNorm;  //!< Time delay RMS fit: Exponential normalization
    double fTDelayRMSExpShift;  //!< Time delay RMS fit: Exponential x shift
    double fTDelayRMSExpScale;  //!< Time delay RMS fit: Exponential scale
    double fNpeScaleNorm;  //!< Npe vs. distance: 1/r^2 scale
    double fNpeScaleShift;  //!< Npe vs. distance: 1/r^2 x shift
    double fQ0;  // Average energy deposited for mips, for charge scaling [GeV]
    double fQPed;  // ADC offset for the single-peak peak mean [ADC]
    double fQSlope;  // Slope in mean ADC / Npe [ADC]
    double fQRMS;  // ADC single-pe spectrum width [ADC]
    double fQThreshold;  //!< ADC charge threshold [ADC]
    double fTResInterpolator;  // Interpolator time resolution [ns]
    double fPropDelay;  // Delay in pulse arrival time [ns/m]
    double fPropDelayError;  // Delay in pulse arrival time, uncertainty [ns/m]
    bool fUseEdep;  //!< Use the true G4 energy deposited, assume mip if false.
    bool fModelTransAtten;  //!< simplify electronics response
    bool fModelLongAtten;  //!< simplify electronics response
    double fAbsLenEff;  // Effective abs. length for transverse Npe scaling [cm]
    double fStripCoincidenceWindow;
    double fTaggerPlaneCoincidenceWindow;
    double fSipmTimeResponse;
    double fCRTClockFreq;
    bool fverbose;
    bool fSumThresh;
    CLHEP::HepRandomEngine& fEngine;

    /**
     * Get the channel trigger time relative to the start of the MC event.
     *
     * @param engine The random number generator engine
     * @param clock The clock to count ticks on
     * @param t0 The starting time (which delay is added to)
     * @param npe Number of observed photoelectrons
     * @param r Distance between the energy deposit and strip readout end [mm]
     * @return The channel trigger time [ns]
*/
    double getChannelTriggerTicks(detinfo::ElecClock& clock,
                                  float t0, float npeMean, float r);

  };

  // Implementation below
  CRTDetSim::CRTDetSim(const fhicl::ParameterSet& pSet)
    : EDProducer{pSet}
    , fG4ModuleLabel{pSet.get<std::string>("G4ModuleLabel","g4")}
    , fGlobalT0Offset{pSet.get<double>("GlobalT0Offset",0.0)}
    , fTDelayNorm{pSet.get<double>("TDelayNorm",4126.)}
    , fTDelayShift{pSet.get<double>("TDelayShift",-300.)}
    , fTDelaySigma{pSet.get<double>("TDelaySigma",90.)}
    , fTDelayOffset{pSet.get<double>("TDelayOffset",-1.5)}
    , fTDelayRMSGausNorm{pSet.get<double>("TDelayRMSGausNorm",2.09)}
    , fTDelayRMSGausShift{pSet.get<double>("TDelayRMSGausShift",7.24)}
    , fTDelayRMSGausSigma{pSet.get<double>("TDelayRMSGausSigma",170.)}
    , fTDelayRMSExpNorm{pSet.get<double>("TDelayRMSExpNorm",1.65)}
    , fTDelayRMSExpShift{pSet.get<double>("TDelayRMSExpShift",75.6)}
    , fTDelayRMSExpScale{pSet.get<double>("TDelayRMSExpScale",79.35)}
    , fNpeScaleNorm{pSet.get<double>("NpeScaleNorm",5261000.)}
    , fNpeScaleShift{pSet.get<double>("NpeScaleShift",-1085.)}
    , fQ0{pSet.get<double>("Q0",0.00175)}
    , fQPed{pSet.get<double>("QPed",63.6)}
    , fQSlope{pSet.get<double>("QSlope",131.9)}
    , fQRMS{pSet.get<double>("QRMS",15.)}
    , fQThreshold{pSet.get<double>("QThreshold",100.0)}
    , fTResInterpolator{pSet.get<double>("TResInterpolator",1.268)}
    , fPropDelay{pSet.get<double>("PropDelay",0.0061)}
    , fPropDelayError{pSet.get<double>("PropDelayError",0.007)}
    , fUseEdep{pSet.get<bool>("UseEdep",true)}
    , fModelTransAtten{pSet.get<bool>("ModelTransAtten",true)}
    , fModelLongAtten{pSet.get<bool>("ModelLongAtten",true)}
    , fAbsLenEff{pSet.get<double>("AbsLenEff",8.5)}
    , fStripCoincidenceWindow{pSet.get<double>("StripCoincidenceWindow", 30.)}
    , fTaggerPlaneCoincidenceWindow{pSet.get<double>("TaggerPlaneCoincidenceWindow", 100.)}
    , fSipmTimeResponse{pSet.get<double>("SipmTimeResponse",2.)}
    , fCRTClockFreq{pSet.get<double>("CRTClockFreq",1.0)}
    , fverbose{pSet.get<bool>("verbose",0)}
    , fSumThresh{pSet.get<bool>("SumThresh",false)}
    , fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "crt", pSet, "Seed"))
  {
    if (fQThreshold<fQPed) { 
      std::cout << " Threshold cannot be lower than pedestal value.  Setting threshold at pedestal" << std::endl;
      fQThreshold=fQPed;
    }

    produces< std::vector<CRTSimData> >();
  }
      
  double CRTDetSim::getChannelTriggerTicks(detinfo::ElecClock& clock,
                                           float t0, float npeMean, float r)
  {
    // Hit timing, with smearing and NPE dependence
    double tDelayMean =
      fTDelayNorm *
        exp(-0.5 * pow((npeMean - fTDelayShift) / fTDelaySigma, 2)) +
      fTDelayOffset;

    double tDelayRMS =
      fTDelayRMSGausNorm *
        exp(-pow(npeMean - fTDelayRMSGausShift, 2) / fTDelayRMSGausSigma) +
      fTDelayRMSExpNorm *
        exp(-(npeMean - fTDelayRMSExpShift) / fTDelayRMSExpScale);

    double tDelay = CLHEP::RandGauss::shoot(&fEngine, tDelayMean, tDelayRMS);

    // Time resolution of the interpolator
    tDelay += CLHEP::RandGauss::shoot(&fEngine, 0, fTResInterpolator);

    // Propagation time
    double tProp = CLHEP::RandGauss::shoot(fPropDelay, fPropDelayError) * r;
    double t = t0 + tProp + tDelay;

    // Get clock ticks
    clock = clock.WithTime(t/1000.);  // SetTime takes microseconds
    // convert to CRT clock ticks instead (1 tick = 1 ns)    
    double cticks = uint(t/fCRTClockFreq);
    

    mf::LogInfo("CRT")
    << "CRT TIMING: t0=" << t0
    << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS
    << ", tDelay=" << tDelay << ", tDelay(interp)="
    << tDelay << ", tProp=" << tProp << ", t=" << t << ", ticks=" << cticks << "\n";

    return cticks;
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

    // not used
    art::ServiceHandle<detinfo::DetectorClocksService> detClocks;
    detinfo::ElecClock trigClock = detClocks->DataFor(evt).TriggerClock();

    // Handle for (truth) AuxDetSimChannels
    art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
    evt.getByLabel(fG4ModuleLabel, channels);

    //access geometry
    geo::GeometryCore const* fGeometryService;            
    fGeometryService = lar::providerFrom<geo::Geometry>();


  ///< pointer to Geometry provider   
    //    art::ServiceHandle<geo::AuxDetGeometry> geoService;
    //    const geo::AuxDetGeometry* geometry = &*geoService;
    // const geo::AuxDetGeometryCore* geoServiceProvider = geometry->GetProviderPtr();

    // Loop through truth AD channels
    for (auto& adsc : *channels) {
      //      const geo::AuxDetGeo& adGeo = geoServiceProvider->AuxDet(adsc.AuxDetID());
      const geo::AuxDetGeo& adGeo = fGeometryService->AuxDet(adsc.AuxDetID());
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
        double eDep = ide.energyDeposited;
	if (eDep>0.0000001) {  // require ~0.01 pe or more to proceed
	double maxEdep = eDep;	
   	int trackID = ide.trackID;
           // Get the hit position in strip's local coordinates
        double x = (ide.entryX + ide.exitX) / 2;
        double y = (ide.entryY + ide.exitY) / 2;
        double z = (ide.entryZ + ide.exitZ) / 2;

        double tTrue = (ide.entryT + ide.exitT) / 2 + fGlobalT0Offset;


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

	//	const geo::AuxDetGeo& adGeo = fGeometryService->AuxDet(adsc.AuxDetID());
	std::string name = adGeo.TotalVolume()->GetName();
	int module,strip;
	sscanf(name.c_str(),"volAuxDet_Module_%d_strip_%d",&module,&strip);
	if (fverbose) std::cout << " channel " << adsc.AuxDetID() << " strip " << strip << 
			" module " << module << " name " << name << std::endl;
	double world[3] = {x, y, z};
	double svHitPosLocal[3];
	adsGeo.WorldToLocal(world, svHitPosLocal);	

        // Distance to the readout end 
	double distToReadout=0;
	if (fModelLongAtten) distToReadout= abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);

	// The expected number of PE, using a quadratic model for the distance
	// dependence, and scaling linearly with deposited energy.
	//  UseEdep flag let's us do a binary response
	double qr = fUseEdep ? 1.0 * eDep / fQ0 : 1.0;
	double npeExpected = (fNpeScaleNorm / pow(distToReadout - fNpeScaleShift, 2) * qr);
	//	if (fverbose) std::cout << " eDep " << eDep << " eDep/fQ0=" << eDep/fQ0 << " npe exp =" << npeExpected << std::endl;

	// Put PE on channels weighted by distance
	double d0 = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);  // L
	double d1 = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]);  // R
	short q0,q1;
	double npeExp0=0;
	double npeExp1=0;
	long npe0 =0; 
	long npe1 =0;
        // Time relative to PPS: set to zero instead of random
        //        uint32_t ppsTicks = CLHEP::RandFlat::shootInt(fEngine, trigClock.Frequency() * 1e6);
        uint32_t ppsTicks = 0;
		
	if (fModelTransAtten) {
	  double abs0 = exp(-d0 / fAbsLenEff);	
	  double abs1 = exp(-d1 / fAbsLenEff);
	  npeExp0 = npeExpected * abs0 / (abs0 + abs1);
	  npeExp1 = npeExpected * abs1 / (abs0 + abs1);
	}
	else {
	  npeExp0 = npeExpected*d0/(d1+d0);
	  npeExp1 = npeExpected*d1/(d1+d0);
	}
	if (fModelTransAtten) {  // should be a different flag to turn on/off smearing
	// Observed PE
          npe0 = CLHEP::RandPoisson::shoot(&fEngine, npeExp0);
          npe1 = CLHEP::RandPoisson::shoot(&fEngine, npeExp1);
	  // SiPM and ADC response: Npe to ADC counts
          q0 = CLHEP::RandGauss::shoot(&fEngine, fQPed + fQSlope * npe0, fQRMS * sqrt(npe0));
          q1 = CLHEP::RandGauss::shoot(&fEngine, fQPed + fQSlope * npe1, fQRMS * sqrt(npe1));
	}
	else {
	  npe0=npeExp0;
	  npe1=npeExp1;
	  q0=fQPed + fQSlope * npe0;
	  q1=fQPed + fQSlope * npe1;
	}

	if (fverbose) std::cout << "npe exp = " << npeExpected << " q0,q1 " << q0 << " " << q1 << 
	  " d0,d1 " << d0 << " " << d1 << std::endl;
      
	// NOT The time relative to trigger in trigger ticks 
	//  trigClock not used
	//  Instead time relative to event time=0 (will include time offset of particles not 
	//       generated at time=0) in CRT ticks (=1 ns)
        uint32_t t0 = getChannelTriggerTicks(trigClock, tTrue, npe0, distToReadout);
        uint32_t t1 = getChannelTriggerTicks(trigClock, tTrue, npe1, distToReadout);

        uint32_t channel0ID = adsc.AuxDetID()*2;
        uint32_t channel1ID = adsc.AuxDetID()*2+1;
	
      // Apply ADC threshold and strip-level coincidence (both fibers fire)
	if (util::absDiff(t0, t1) < fStripCoincidenceWindow) {
	  if ( (q0 > fQThreshold && q1 > fQThreshold) || ((fSumThresh) && (q0+q1)>fQThreshold)) {
	    crtHits->push_back(CRTSimData(channel0ID, t0, ppsTicks, q0, trackID));
	    crtHits->push_back(CRTSimData(channel1ID, t1, ppsTicks, q1, trackID));
	    if (fverbose) std::cout << "both written to event q0: " << q0 << " q1 " << q1 <<std::endl;
	  }
	  else {
	    if (fverbose) std::cout << "below threshold: q0 "  << q0 << " q1 " << q1 <<std::endl;
	  }
	} // if coincidence
	else {
	  if (fverbose) std::cout << "failed time coincidence requirement: t0  " << t0 << " t1 " << t1 << std::endl;
	}  // else coincidence
	} //zero energy threshold
    } //loop over ides
    }// loop over channels

    evt.put(std::move(crtHits));
  }
  

  DEFINE_ART_MODULE(CRTDetSim)
}  //namespace CRT
