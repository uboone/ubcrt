////////////////////////////////////////////////////////////////////////
// Class:       G4VetoSignalImpact
// Plugin Type: analyzer (art v2_11_03)
// File:        G4VetoSignalImpact_module.cc
//
// Generated at Wed Oct  3 14:42:27 2018 by David Caratelli using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

// Backtrack which particles come from the neutrino
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
//#include "larsim/MCCheater/BackTrackerService.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

//#include "LLBasicTool/GeoAlgo/GeoAlgo.h"

class G4VetoSignalImpact;

class G4VetoSignalImpact : public art::EDAnalyzer
{
public:
  explicit G4VetoSignalImpact(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  G4VetoSignalImpact(G4VetoSignalImpact const &) = delete;
  G4VetoSignalImpact(G4VetoSignalImpact &&) = delete;
  G4VetoSignalImpact &operator=(G4VetoSignalImpact const &) = delete;
  G4VetoSignalImpact &operator=(G4VetoSignalImpact &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  //void endSubRun(const art::SubRun &sr);

  art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const &e, std::string _geant_producer, int geant_track_id);

private:
  // Declare member data here.

  std::pair<bool, TVector3> CRTIntersection(const simb::MCTrajectory &traj);
  std::pair<bool, TVector3> InCRT(const TLorentzVector &p1, const TLorentzVector &p2);

  // set CRT geometry configuration here
  double fThickness;                         // CRT panel thickness (assumed the same for all panels...)
  double fTx1, fTx2, fTy1, fTy2, fTz1, fTz2; // assuming each panel is a rectangle. 6 coordinates is all that is needed to define the geometry.
  double fBx1, fBx2, fBy1, fBy2, fBz1, fBz2; // assuming each panel is a rectangle. 6 coordinates is all that is needed to define the geometry.
  double fAx1, fAx2, fAy1, fAy2, fAz1, fAz2; // assuming each panel is a rectangle. 6 coordinates is all that is needed to define the geometry.
  double fCx1, fCx2, fCy1, fCy2, fCz1, fCz2; // assuming each panel is a rectangle. 6 coordinates is all that is needed to define the geometry.

  TTree *_pottree;
  float _pot;
  int _run_sr;
  float _subrun_sr;
  // variables for TTree
  TTree *_tree;
  int _intersections;
  std::vector<float> _x_cross;
  std::vector<float> _y_cross;
  std::vector<float> _z_cross;
  // information on particle crossing the CRT:
  // PDG code, energy @ production, start and end point
  std::vector<int> _cross_pdg;
  std::vector<float> _cross_momentum;
  std::vector<float> _cross_xstart;
  std::vector<float> _cross_ystart;
  std::vector<float> _cross_zstart;
  std::vector<float> _cross_xend;
  std::vector<float> _cross_yend;
  std::vector<float> _cross_zend;
  // neutrino fields if there is a neutrino:
  uint fNum_nu;
  std::vector<float> fNu_vtx_x;
  std::vector<float> fNu_vtx_y;
  std::vector<float> fNu_vtx_z;
  std::vector<float> fNu_E;
  std::vector<float> fNu_time;
  std::vector<int> fNu_pdg_code;
  std::vector<bool> fNu_ccnc;
};

G4VetoSignalImpact::G4VetoSignalImpact(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{

  fThickness = p.get<double>("Thickness");
  fTx1 = p.get<double>("Tx1");
  fTy1 = p.get<double>("Ty1");
  fTz1 = p.get<double>("Tz1");
  fTx2 = p.get<double>("Tx2");
  fTy2 = p.get<double>("Ty2");
  fTz2 = p.get<double>("Tz2");
  fBx1 = p.get<double>("Bx1");
  fBy1 = p.get<double>("By1");
  fBz1 = p.get<double>("Bz1");
  fBx2 = p.get<double>("Bx2");
  fBy2 = p.get<double>("By2");
  fBz2 = p.get<double>("Bz2");
  fAx1 = p.get<double>("Ax1");
  fAy1 = p.get<double>("Ay1");
  fAz1 = p.get<double>("Az1");
  fAx2 = p.get<double>("Ax2");
  fAy2 = p.get<double>("Ay2");
  fAz2 = p.get<double>("Az2");
  fCx1 = p.get<double>("Cx1");
  fCy1 = p.get<double>("Cy1");
  fCz1 = p.get<double>("Cz1");
  fCx2 = p.get<double>("Cx2");
  fCy2 = p.get<double>("Cy2");
  fCz2 = p.get<double>("Cz2");
}

void G4VetoSignalImpact::analyze(art::Event const &e)
{

  // this module assumes all MCParticles are from neutrino interactions (i.e. genie simulation only, no corsika)

  // the goal of the module is to establish if any charged partciles deposit energy in a CRT module and save details about these interactions.
  // information to be stored:
  // neutrino energy, interaction type (CC, NC), interaction position (in TPC, in cryo, dirt)
  // details of particle crossing CRT (e-, muon, etc...)
  // CRT panel(s) hit

  _intersections = 0;
  fNum_nu = 0;
  fNu_vtx_x.clear();
  fNu_vtx_y.clear();
  fNu_vtx_z.clear();
  fNu_E.clear();
  fNu_time.clear();
  fNu_pdg_code.clear();
  fNu_ccnc.clear();

  _x_cross.clear();
  _y_cross.clear();
  _z_cross.clear();

  _cross_pdg.clear();
  _cross_momentum.clear();
  _cross_xstart.clear();
  _cross_ystart.clear();
  _cross_zstart.clear();
  _cross_xend.clear();
  _cross_yend.clear();
  _cross_zend.clear();

  // load neutrino mctruth
  auto const &generator_handle = e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  auto const &generator(*generator_handle);
  fNum_nu = generator.size();

  for (auto &gen : generator)
  {

    if (gen.Origin() == simb::kBeamNeutrino)
    {
      fNu_pdg_code.push_back(gen.GetNeutrino().Nu().PdgCode());
      fNu_time.push_back(gen.GetNeutrino().Nu().T());
      fNu_E.push_back(gen.GetNeutrino().Nu().E());
      fNu_ccnc.push_back(gen.GetNeutrino().CCNC());
      fNu_vtx_x.push_back(gen.GetNeutrino().Nu().Vx());
      fNu_vtx_y.push_back(gen.GetNeutrino().Nu().Vy());
      fNu_vtx_z.push_back(gen.GetNeutrino().Nu().Vz());
    }
  }

  // load MCParticles
  art::Handle<std::vector<simb::MCParticle>> MCParticle_h;
  e.getByLabel("largeant", MCParticle_h);

  if (!MCParticle_h.isValid())
  {
    std::cerr << "\033[93m[ERROR]\033[00m ... could not locate MC Particles!" << std::endl;
    throw std::exception();
  }

  for (size_t n = 0; n < MCParticle_h->size(); n++)
  {

    auto const &mcpart = MCParticle_h->at(n);

    // is this particle neutral? if so, ignore
    // specifically, exclude photons and neutrons)
    uint pdg = abs(mcpart.PdgCode());
    bool pdg_ok = pdg == 11 or pdg == 13 or pdg == 211 or pdg == 111 or pdg == 2212;
    if (!pdg_ok or mcpart.E()<0.05)
      continue;

    // is this particle related to the neutrino interaction?
    const art::Ptr<simb::MCTruth> mctruth = TrackIDToMCTruth(e, "largeant", mcpart.TrackId());
    if (mctruth->Origin() != simb::kBeamNeutrino)
      continue;

    // ask for whether this particle intersects the CRT
    auto crtintersection = CRTIntersection(mcpart.Trajectory());

    if (crtintersection.first == true)
    {

      _intersections += 1;
      _x_cross.push_back(crtintersection.second.X());
      _y_cross.push_back(crtintersection.second.Y());
      _z_cross.push_back(crtintersection.second.Z());

      _cross_pdg.push_back(mcpart.PdgCode());
      int nstep = mcpart.Trajectory().size() - 1;
      _cross_momentum.push_back(mcpart.Trajectory().Momentum(0).Vect().Mag());
      _cross_xstart.push_back(mcpart.Trajectory().Position(0).X());
      _cross_ystart.push_back(mcpart.Trajectory().Position(0).Y());
      _cross_zstart.push_back(mcpart.Trajectory().Position(0).Z());
      _cross_xend.push_back(mcpart.Trajectory().Position(nstep).X());
      _cross_yend.push_back(mcpart.Trajectory().Position(nstep).Y());
      _cross_zend.push_back(mcpart.Trajectory().Position(nstep).Z());
    }

  } // for all MCParticles in the event

  std::cout << "[G4VetoSignalImpact] True neutrinos found: " << fNum_nu << "\tIntersections found: " << _intersections << std::endl;
  _tree->Fill();

  return;
}

std::pair<bool, TVector3> G4VetoSignalImpact::CRTIntersection(const simb::MCTrajectory &traj)
{

  for (size_t t = 1; t < traj.size(); t++)
  {

    auto intersection = InCRT(traj.Position(t - 1), traj.Position(t));

    if (intersection.first == false)
      continue;

    // if we get to this point -> the particle intersects the CRT
    return std::make_pair(true, intersection.second);

  } // for all steps along the trajectory

  return std::make_pair(false, TVector3(0., 0., 0.));
}

std::pair<bool, TVector3> G4VetoSignalImpact::InCRT(const TLorentzVector &p1, const TLorentzVector &p2)
{

  // get points associated with the TLorentzVectors:

  auto pt1 = p1.Vect();
  auto pt2 = p2.Vect();

  // check if the point intersects the top panel
  double t = (fTy1 - pt1.Y()) / (pt2.Y() - pt1.Y());
  // if t < 0 or > 1 then the intersection is beyond the segment
  if ((t > 0) && (t <= 1))
  {
    double ptTx = pt1.X() + (pt2.X() - pt1.X()) * t;
    double ptTz = pt1.Z() + (pt2.Z() - pt1.Z()) * t;
    if ((ptTx > fTx1) && (ptTx < fTx2) && (ptTz > fTz1) && (ptTz < fTz2))
    {
      double ptTy = pt1.Y() + (pt2.Y() - pt1.Y()) * t;
      return std::make_pair(true, TVector3(ptTx, ptTy, ptTz));
    } // if they intersec
  }   // if t is between 0 and 1

  // check if the point intersects the bottom panel
  t = (fBy1 - pt1.Y()) / (pt2.Y() - pt1.Y());
  // if t < 0 or > 1 then the intersection is beyond the segment
  if ((t > 0) && (t <= 1))
  {
    double ptBx = pt1.X() + (pt2.X() - pt1.X()) * t;
    double ptBz = pt1.Z() + (pt2.Z() - pt1.Z()) * t;
    if ((ptBx > fBx1) && (ptBx < fBx2) && (ptBz > fBz1) && (ptBz < fBz2))
    {
      double ptBy = pt1.Y() + (pt2.Y() - pt1.Y()) * t;
      return std::make_pair(true, TVector3(ptBx, ptBy, ptBz));
    } // if they intersec
  }   // if t is between 0 and 1

  // check if the point intersects the anode panel
  t = (fAy1 - pt1.Y()) / (pt2.Y() - pt1.Y());
  // if t < 0 or > 1 then the intersection is beyond the segment
  if ((t > 0) && (t <= 1))
  {
    double ptAx = pt1.X() + (pt2.X() - pt1.X()) * t;
    double ptAz = pt1.Z() + (pt2.Z() - pt1.Z()) * t;
    if ((ptAx > fAx1) && (ptAx < fAx2) && (ptAz > fAz1) && (ptAz < fAz2))
    {
      double ptAy = pt1.Y() + (pt2.Y() - pt1.Y()) * t;
      return std::make_pair(true, TVector3(ptAx, ptAy, ptAz));
    } // if they intersec
  }   // if t is between 0 and 1

  // check if the point intersects the cathode panel
  t = (fCy1 - pt1.Y()) / (pt2.Y() - pt1.Y());
  // if t < 0 or > 1 then the intersection is beyond the segment
  if ((t > 0) && (t <= 1))
  {
    double ptCx = pt1.X() + (pt2.X() - pt1.X()) * t;
    double ptCz = pt1.Z() + (pt2.Z() - pt1.Z()) * t;
    if ((ptCx > fCx1) && (ptCx < fCx2) && (ptCz > fCz1) && (ptCz < fCz2))
    {
      double ptCy = pt1.Y() + (pt2.Y() - pt1.Y()) * t;
      return std::make_pair(true, TVector3(ptCx, ptCy, ptCz));
    } // if they intersec
  }   // if t is between 0 and 1

  return std::make_pair(false, TVector3(0., 0., 0.));
}

void G4VetoSignalImpact::beginJob()
{
  // Implementation of optional member function here.

  art::ServiceHandle<art::TFileService> tfs;
  //// Tree for every subrun
  _pottree = tfs->make<TTree>("pot", "POT Tree");
  _pottree->Branch("run", &_run_sr, "run/i");
  _pottree->Branch("subrun", &_subrun_sr, "subrun/i");
  _pottree->Branch("pot", &_pot, "pot/F");

  _tree = tfs->make<TTree>("_tree", "G4 CRT Veto");
  _tree->Branch("_intersections", &_intersections, "intersections/I");
  _tree->Branch("num_nu", &fNum_nu, "num_nu/I");
  _tree->Branch("nu_vtx_x", "std::vector< float >", &fNu_vtx_x);
  _tree->Branch("nu_vtx_y", "std::vector< float >", &fNu_vtx_y);
  _tree->Branch("nu_vtx_z", "std::vector< float >", &fNu_vtx_z);
  _tree->Branch("nu_E", "std::vector< float >", &fNu_E);
  _tree->Branch("nu_time", "std::vector< float >", &fNu_time);
  _tree->Branch("nu_pdg_code", "std::vector< int >", &fNu_pdg_code);
  _tree->Branch("nu_ccnc", "std::vector< bool >", &fNu_ccnc);

  _tree->Branch("_x_cross", "std::vector< float >", &_x_cross);
  _tree->Branch("_y_cross", "std::vector< float >", &_y_cross);
  _tree->Branch("_z_cross", "std::vector< float >", &_z_cross);
  _tree->Branch("_cross_pdg", "std::vector< int >", &_cross_pdg);
  _tree->Branch("_cross_momentum", "std::vector< float >", &_cross_momentum);
  _tree->Branch("_cross_xstart", "std::vector< float >", &_cross_xstart);
  _tree->Branch("_cross_ystart", "std::vector< float >", &_cross_ystart);
  _tree->Branch("_cross_zstart", "std::vector< float >", &_cross_zstart);
  _tree->Branch("_cross_xend", "std::vector< float >", &_cross_xend);
  _tree->Branch("_cross_yend", "std::vector< float >", &_cross_yend);
  _tree->Branch("_cross_zend", "std::vector< float >", &_cross_zend);
}

void G4VetoSignalImpact::endJob()
{
  // Implementation of optional member function here.
}

art::Ptr<simb::MCTruth> G4VetoSignalImpact::TrackIDToMCTruth(art::Event const &e, std::string _geant_producer, int geant_track_id)
{

  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;

  lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, truthToParticles, particlesToTruth);

  for (auto iter : particlesToTruth)
  {
    if (iter.first->TrackId() == geant_track_id)
    {
      return iter.second;
    }
  }

  art::Ptr<simb::MCTruth> null_ptr;
  return null_ptr;
}

/*
void G4VetoSignalImpact::endSubRun(const art::SubRun &sr)
{

  _run_sr = sr.run();
  _subrun_sr = sr.subRun();

  art::Handle<sumdata::POTSummary> potListHandle;
  if (!m_is_data)
  {
    if (sr.getByLabel("generator", potListHandle))
    {
      _pot = potListHandle->totpot;
      std::cout << "[Prototyping::endSubRun] POT for SubRun: " << fPot << std::endl;
    }
    else
      _pot = 0.;
  }
  else
  {
    if (sr.getByLabel("beamdata", "bnbETOR860", potListHandle))
      _pot = potListHandle->totpot;
    else
      _pot = 0.;
  }

  _pottree->Fill();

  if (m_verb)
  {
    std::cout << "string_process has mamebers: " << string_process.size() << std::endl;
    for (auto elem : string_process)
    {
      std::cout << elem << ", ";
    }
  }
}
*/

DEFINE_ART_MODULE(G4VetoSignalImpact)
