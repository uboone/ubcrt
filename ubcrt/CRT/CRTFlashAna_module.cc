////////////////////////////////////////////////////////////////////////
// Class:       CRTFlashAna
// Plugin Type: analyzer (art v2_05_01)
// File:        CRTFlashAna_module.cc
//
// Generated at Fri Jul  6 13:52:52 2018 by Herbert Greenlee using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"
#include "uboone/CRT/CRTProducts/CRTHit.hh"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1F.h"

namespace crt {
  class CRTFlashAna;
}


class crt::CRTFlashAna : public art::EDAnalyzer {
public:
  explicit CRTFlashAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTFlashAna(CRTFlashAna const &) = delete;
  CRTFlashAna(CRTFlashAna &&) = delete;
  CRTFlashAna & operator = (CRTFlashAna const &) = delete;
  CRTFlashAna & operator = (CRTFlashAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Methods.

  void add_algorithm(const std::string& algo);

private:

  // Declare member data here.

  // Fcl parameters.

  std::string fDAQTimeModuleLabel;
  std::string fTriggerModuleLabel;
  std::string fSWTriggerModuleLabel;
  std::string fFlashModuleLabel;
  std::string fCRTModuleLabel;
  double fFlashMinPE;

  // Directories.

  art::TFileDirectory fTopDir;
  std::set<std::string> fAlgos;

  // Histograms.

  std::map<std::string, TH1F*> fHPE;        // Flash PE.
  std::map<std::string, TH1F*> fHcrt0;      // Flash vs. CRT t0 time difference.
  std::map<std::string, TH1F*> fHcrt0x;     // Flash vs. CRT t0 time difference expanded.
  std::map<std::string, TH1F*> fHcrt0d;     // Flash vs. CRT t0 time difference detail.
  std::map<std::string, TH1F*> fHcrtadj0;   // Flash vs. CRT t0 (adj) time difference.
  std::map<std::string, TH1F*> fHcrtadj0x;  // Flash vs. CRT t0 (adj) time difference expanded.
  std::map<std::string, TH1F*> fHcrtadj0d;  // Flash vs. CRT t0 (adj) time difference detail.
};

// Constructor.

crt::CRTFlashAna::CRTFlashAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fDAQTimeModuleLabel(p.get<std::string>("DAQTimeModuleLabel")),
  fTriggerModuleLabel(p.get<std::string>("TriggerModuleLabel")),
  fSWTriggerModuleLabel(p.get<std::string>("SWTriggerModuleLabel")),
  fFlashModuleLabel(p.get<std::string>("FlashModuleLabel")),
  fCRTModuleLabel(p.get<std::string>("CRTModuleLabel")),
  fFlashMinPE(p.get<double>("FlashMinPE")),
  fTopDir(art::ServiceHandle<art::TFileService>()->mkdir("CRTFlashAna", "CRT vs. Flash Timing"))
{
  std::string all("All");
  add_algorithm(all);
}

void crt::CRTFlashAna::analyze(art::Event const & evt)
{
  // Get GPS time.

  art::Handle<raw::DAQHeaderTimeUBooNE> htime;
  evt.getByLabel(fDAQTimeModuleLabel, htime);
  if(!htime.isValid()) {
    std::cout << "No DAQHeaderTimeUBooNE data product with module label " 
	      << fDAQTimeModuleLabel << std::endl;
    return;
  }
  time_t gps_time = htime->gps_time();
  unsigned int gps_sec = (gps_time >> 32);
  unsigned int gps_nsec = (gps_time & 0xffffffff);

  time_t gps_adj_time = htime->gps_adj_time();
  unsigned int gps_adj_sec = (gps_adj_time >> 32);
  unsigned int gps_adj_nsec = (gps_adj_time & 0xffffffff);

  time_t ntp_time = htime->ntp_time();
  unsigned int ntp_sec = (ntp_time >> 32);
  unsigned int ntp_nsec = (ntp_time & 0xffffffff);

  std::cout << "GPS time =     "
	    << gps_sec << " seconds, "
	    << gps_nsec << " nanoseconds." << std::endl;
  std::cout << "GPS adj time = "
	    << gps_adj_sec << " seconds, "
	    << gps_adj_nsec << " nanoseconds." << std::endl;
  std::cout << "NTP time =     "
	    << ntp_sec << " seconds, "
	    << ntp_nsec << " nanoseconds." << std::endl;

  // Get hardware trigger information.

  art::Handle<std::vector<raw::Trigger> > htrig;
  evt.getByLabel(fTriggerModuleLabel, htrig);
  if(!htrig.isValid()) {
    std::cout << "No Trigger data product with module label " 
	      << fTriggerModuleLabel << std::endl;
    return;
  }
  for(auto const& trig : *htrig) {
    double trigger_time = trig.TriggerTime();
    std::cout << "Trigger time = " << trigger_time << std::endl;
  }

  // Get Flashes.

  art::Handle<std::vector<recob::OpFlash> > hflash;
  evt.getByLabel(fFlashModuleLabel, hflash);
  if(!hflash.isValid())
    return;

  // Get CRT hits.

  art::Handle<std::vector<crt::CRTHit> > hcrthit;
  evt.getByLabel(fCRTModuleLabel, hcrthit);
  if(!hcrthit.isValid())
    return;

  // Get software trigger information.

  std::string all("All");
  art::Handle<raw::ubdaqSoftwareTriggerData> hswtrig;
  evt.getByLabel(fSWTriggerModuleLabel, hswtrig);
  if(!hswtrig.isValid()) {
    std::cout << "No ubdaqSoftwareTriggerData data product with module label " 
	      << fSWTriggerModuleLabel << std::endl;
    return;
  }

  // Loop over software triggers.

  auto const& swtrig = *hswtrig; 
  std::vector<std::string> algos = swtrig.getListOfAlgorithms();
  for(auto const& algo : algos) {
    bool pass_algo = swtrig.passedAlgo(algo);
    bool pass_prescale = swtrig.passedPrescaleAlgo(algo);
    std::cout << algo << " " << pass_algo << " " << pass_prescale << std::endl;
    if(pass_algo && pass_prescale) {
      add_algorithm(algo);

      // Make a double loop over CRT hits and flashes.

      for(auto const& flash : *hflash) {
	double flash_time = flash.Time();   // microseconds.
	double flash_pe = flash.TotalPE();  // Total PE.
	fHPE[algo]->Fill(flash_pe);

	if(flash_pe >= fFlashMinPE) {
  
	  for(auto const& crthit : *hcrthit) {

	    // Calculate CRT hit time relative to gps time (units microseconds).

	    double crthit_time =
	      1.e6 * (double(crthit.ts0_s) - double(gps_sec)) +
	      1.e-3 * (double(crthit.ts0_ns) - double(gps_nsec));
	    double crthit_adj_time =
	      1.e6 * (double(crthit.ts0_s) - double(gps_adj_sec)) +
	      1.e-3 * (double(crthit.ts0_ns) - double(gps_adj_nsec));

	    // Fill histograms.

	    fHcrt0[algo]->Fill(crthit_time - flash_time);
	    fHcrt0x[algo]->Fill(crthit_time - flash_time);
	    fHcrt0d[algo]->Fill(crthit_time - flash_time);
	    fHcrt0[all]->Fill(crthit_time - flash_time);
	    fHcrt0x[all]->Fill(crthit_time - flash_time);
	    fHcrt0d[all]->Fill(crthit_time - flash_time);
	    fHcrtadj0[algo]->Fill(crthit_adj_time - flash_time);
	    fHcrtadj0x[algo]->Fill(crthit_adj_time - flash_time);
	    fHcrtadj0d[algo]->Fill(crthit_adj_time - flash_time);
	    fHcrtadj0[all]->Fill(crthit_adj_time - flash_time);
	    fHcrtadj0x[all]->Fill(crthit_adj_time - flash_time);
	    fHcrtadj0d[all]->Fill(crthit_adj_time - flash_time);
	  }
	}
      }
    }
  }
}

void crt::CRTFlashAna::add_algorithm(const std::string& algo)
{
  if(fAlgos.count(algo) == 0) {
    fAlgos.insert(algo);
    art::TFileDirectory dir = fTopDir.mkdir(algo);
    fHPE[algo] = dir.make<TH1F>("FlashPE", "Flash PE", 100, 0., 1000.);
    fHPE[algo]->GetXaxis()->SetTitle("Flash PE (ADC)");
    fHcrt0[algo] = dir.make<TH1F>("crt0", "CRT vs. Flash Time Difference", 1000, -5000., 5000.);
    fHcrt0[algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
    fHcrt0x[algo] = dir.make<TH1F>("crt0x", "CRT vs. Flash Time Difference Expanded",
				  500, -500., 0.);
    fHcrt0x[algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
    fHcrt0d[algo] = dir.make<TH1F>("crt0d", "CRT vs. Flash Time Difference Detail",
				  100, -100., -50.);
    fHcrt0d[algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
    fHcrtadj0[algo] = dir.make<TH1F>("crtadj0", "CRT vs. Flash Time Difference (Adj)",
				     1000, -5000., 5000.);
    fHcrtadj0[algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
    fHcrtadj0x[algo] = dir.make<TH1F>("crtadj0x", "CRT vs. Flash Time Difference (Adj) Expanded",
				     500, -500., 0.);
    fHcrtadj0x[algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
    fHcrtadj0d[algo] = dir.make<TH1F>("crtadj0d", "CRT vs. Flash Time Difference (Adj) Detail",
				     100, -100., -50.);
    fHcrtadj0d[algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
  }
}

DEFINE_ART_MODULE(crt::CRTFlashAna)
