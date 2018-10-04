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

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

//#include "LLBasicTool/GeoAlgo/GeoAlgo.h"

class G4VetoSignalImpact;


class G4VetoSignalImpact : public art::EDAnalyzer {
public:
  explicit G4VetoSignalImpact(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  G4VetoSignalImpact(G4VetoSignalImpact const &) = delete;
  G4VetoSignalImpact(G4VetoSignalImpact &&) = delete;
  G4VetoSignalImpact & operator = (G4VetoSignalImpact const &) = delete;
  G4VetoSignalImpact & operator = (G4VetoSignalImpact &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  std::pair<bool, TVector3> CRTIntersection(const simb::MCTrajectory& traj);
  std::pair<bool, TVector3> InCRT(const TLorentzVector& p1, const TLorentzVector& p2);

  // set CRT geometry configuration here
  double fThickness; // CRT panel thickness (assumed the same for all panels...)
  double fTx1, fTx2, fTy1, fTy2, fTz1, fTz2; // assuming each panel is a rectangle. 6 coordinates is all that is needed to define the geometry.
  double fBx1, fBx2, fBy1, fBy2, fBz1, fBz2; // assuming each panel is a rectangle. 6 coordinates is all that is needed to define the geometry.
  double fAx1, fAx2, fAy1, fAy2, fAz1, fAz2; // assuming each panel is a rectangle. 6 coordinates is all that is needed to define the geometry.
  double fCx1, fCx2, fCy1, fCy2, fCz1, fCz2; // assuming each panel is a rectangle. 6 coordinates is all that is needed to define the geometry.

  // variables for TTree
  TTree* _tree;
  int _intersections;
  double _nu_e;
  int _nu_ccnc;
  double _nu_vtx_x;
  double _nu_vtx_y;
  double _nu_vtx_z;
  double _x_cross;
  double _y_cross;
  double _z_cross;
  // information on particle crossing the CRT:
  // PDG code, energy @ production, start and end point
  int    _cross_pdg;
  double _cross_momentum;
  double _cross_xstart;
  double _cross_ystart;
  double _cross_zstart;
  double _cross_xend;
  double _cross_yend;
  double _cross_zend;
  
};


G4VetoSignalImpact::G4VetoSignalImpact(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  fThickness = p.get<double>("Thickness");
  fTx1       = p.get<double>("Tx1");
  fTy1       = p.get<double>("Ty1");
  fTz1       = p.get<double>("Tz1");
  fTx2       = p.get<double>("Tx2");
  fTy2       = p.get<double>("Ty2");
  fTz2       = p.get<double>("Tz2");
  fBx1       = p.get<double>("Bx1");
  fBy1       = p.get<double>("By1");
  fBz1       = p.get<double>("Bz1");
  fBx2       = p.get<double>("Bx2");
  fBy2       = p.get<double>("By2");
  fBz2       = p.get<double>("Bz2");
  fAx1       = p.get<double>("Ax1");
  fAy1       = p.get<double>("Ay1");
  fAz1       = p.get<double>("Az1");
  fAx2       = p.get<double>("Ax2");
  fAy2       = p.get<double>("Ay2");
  fAz2       = p.get<double>("Az2");
  fCx1       = p.get<double>("Cx1");
  fCy1       = p.get<double>("Cy1");
  fCz1       = p.get<double>("Cz1");
  fCx2       = p.get<double>("Cx2");
  fCy2       = p.get<double>("Cy2");
  fCz2       = p.get<double>("Cz2");

}

void G4VetoSignalImpact::analyze(art::Event const & e)
{

  // this module assumes all MCParticles are from neutrino interactions (i.e. genie simulation only, no corsika)

  // the goal of the module is to establish if any charged partciles deposit energy in a CRT module and save details about these interactions.
  // information to be stored:
  // neutrino energy, interaction type (CC, NC), interaction position (in TPC, in cryo, dirt)
  // details of particle crossing CRT (e-, muon, etc...)
  // CRT panel(s) hit
  
  _intersections = 0;
  _nu_e          = 0;
  _nu_ccnc       = 0;
  _nu_vtx_x      = 0;
  _nu_vtx_y      = 0;
  _nu_vtx_z      = 0;
  _x_cross       = 0;
  _y_cross       = 0;
  _z_cross       = 0;

  _cross_pdg = 0;
  _cross_momentum = 0;
  _cross_xstart = 0;
  _cross_ystart = 0;
  _cross_zstart = 0;
  _cross_xend   = 0;
  _cross_yend   = 0;
  _cross_zend   = 0;

  // load neutrino mctruth
  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");

  auto mct = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();

  auto nu     = neutrino.Nu();
  auto lepton = neutrino.Lepton();
  auto ccnc   = neutrino.CCNC();

  _nu_e = nu.Trajectory().E(0);
  _nu_ccnc = ccnc;

  // vertex coordinates from neutrino end point
  _nu_vtx_x = nu.Trajectory().X( nu.Trajectory().size() - 1 );
  _nu_vtx_y = nu.Trajectory().Y( nu.Trajectory().size() - 1 );
  _nu_vtx_z = nu.Trajectory().Z( nu.Trajectory().size() - 1 );

  // load MCParticles
  art::Handle< std::vector<simb::MCParticle> > MCParticle_h;
  e.getByLabel( "largeant", MCParticle_h );

  if(!MCParticle_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate MC Particles!"<<std::endl;
    throw std::exception();
  }


  for (size_t n=0; n < MCParticle_h->size(); n++) {

    auto const& mcpart = MCParticle_h->at(n);

    // is this particle neutral? if so, ignore
    // specifically, exclude photons and neutrons)
    if ( (mcpart.PdgCode() == 22) || (mcpart.PdgCode() == 2112) || (fabs(mcpart.PdgCode()) == 14) || (fabs(mcpart.PdgCode()) == 12) )
      continue;
    
    // ask for whether this particle intersects the CRT
    auto crtintersection =  CRTIntersection(mcpart.Trajectory());

    if (crtintersection.first == true) {
      
      _intersections += 1;
      _x_cross = crtintersection.second.X();
      _y_cross = crtintersection.second.Y();
      _z_cross = crtintersection.second.Z();

      _cross_pdg = mcpart.PdgCode();
      int nstep = mcpart.Trajectory().size()-1;
      _cross_momentum = mcpart.Trajectory().Momentum(0).Vect().Mag();
      _cross_xstart = mcpart.Trajectory().Position(0).X();
      _cross_ystart = mcpart.Trajectory().Position(0).Y();
      _cross_zstart = mcpart.Trajectory().Position(0).Z();
      _cross_xend = mcpart.Trajectory().Position(nstep).X();
      _cross_yend = mcpart.Trajectory().Position(nstep).Y();
      _cross_zend = mcpart.Trajectory().Position(nstep).Z();

    }
    
  }// for all MCParticles in the event

  _tree->Fill();
  
  return;
}

std::pair<bool, TVector3> G4VetoSignalImpact::CRTIntersection(const simb::MCTrajectory& traj) {
  
  for (size_t t=1; t < traj.size(); t++) {

    auto intersection = InCRT(traj.Position(t-1),traj.Position(t));

    if (intersection.first == false)
      continue;

    // if we get to this point -> the particle intersects the CRT
    return std::make_pair(true, intersection.second);
    
  }// for all steps along the trajectory

  return std::make_pair(false, TVector3(0.,0.,0.) );
}

std::pair<bool, TVector3> G4VetoSignalImpact::InCRT(const TLorentzVector& p1, const TLorentzVector& p2) {

  // get points associated with the TLorentzVectors:

  auto pt1 = p1.Vect();
  auto pt2 = p2.Vect();

  // check if the point intersects the top panel
  double t = (fTy1 - pt1.Y())/(pt2.Y()-pt1.Y());
  // if t < 0 or > 1 then the intersection is beyond the segment
  if ( (t > 0) && (t <= 1)) {
    double ptTx = pt1.X() + (pt2.X()-pt1.X())*t;
    double ptTz = pt1.Z() + (pt2.Z()-pt1.Z())*t;
    if ( (ptTx > fTx1) && (ptTx < fTx2) && (ptTz > fTz1) && (ptTz < fTz2) ) {
      double ptTy = pt1.Y() + (pt2.Y()-pt1.Y())*t;
      return std::make_pair(true, TVector3(ptTx,ptTy,ptTz) );
    }// if they intersec
  }// if t is between 0 and 1

  // check if the point intersects the bottom panel
  t = (fBy1 - pt1.Y())/(pt2.Y()-pt1.Y());
  // if t < 0 or > 1 then the intersection is beyond the segment
  if ( (t > 0) && (t <= 1)) {
    double ptBx = pt1.X() + (pt2.X()-pt1.X())*t;
    double ptBz = pt1.Z() + (pt2.Z()-pt1.Z())*t;
    if ( (ptBx > fBx1) && (ptBx < fBx2) && (ptBz > fBz1) && (ptBz < fBz2) ) {
      double ptBy = pt1.Y() + (pt2.Y()-pt1.Y())*t;
      return std::make_pair(true, TVector3(ptBx,ptBy,ptBz) );
    }// if they intersec
  }// if t is between 0 and 1

  // check if the point intersects the anode panel
  t = (fAy1 - pt1.Y())/(pt2.Y()-pt1.Y());
  // if t < 0 or > 1 then the intersection is beyond the segment
  if ( (t > 0) && (t <= 1)) {
    double ptAx = pt1.X() + (pt2.X()-pt1.X())*t;
    double ptAz = pt1.Z() + (pt2.Z()-pt1.Z())*t;
    if ( (ptAx > fAx1) && (ptAx < fAx2) && (ptAz > fAz1) && (ptAz < fAz2) ) {
      double ptAy = pt1.Y() + (pt2.Y()-pt1.Y())*t;
      return std::make_pair(true, TVector3(ptAx,ptAy,ptAz) );
    }// if they intersec
  }// if t is between 0 and 1

  // check if the point intersects the cathode panel
  t = (fCy1 - pt1.Y())/(pt2.Y()-pt1.Y());
  // if t < 0 or > 1 then the intersection is beyond the segment
  if ( (t > 0) && (t <= 1)) {
    double ptCx = pt1.X() + (pt2.X()-pt1.X())*t;
    double ptCz = pt1.Z() + (pt2.Z()-pt1.Z())*t;
    if ( (ptCx > fCx1) && (ptCx < fCx2) && (ptCz > fCz1) && (ptCz < fCz2) ) {
      double ptCy = pt1.Y() + (pt2.Y()-pt1.Y())*t;
      return std::make_pair(true, TVector3(ptCx,ptCy,ptCz) );
    }// if they intersec
  }// if t is between 0 and 1

  return std::make_pair(false, TVector3(0.,0.,0.) );
}

void G4VetoSignalImpact::beginJob()
{
  // Implementation of optional member function here.

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","G4 CRT Veto");
  _tree->Branch("_intersections",&_intersections,"intersections/I");
  _tree->Branch("_nu_e",&_nu_e,"nu_e/D");
  _tree->Branch("_nu_ccnc",&_nu_ccnc,"nu_ccnc/I");
  _tree->Branch("_nu_vtx_x",&_nu_vtx_x,"nu_vtx_x/D");
  _tree->Branch("_nu_vtx_y",&_nu_vtx_y,"nu_vtx_y/D");
  _tree->Branch("_nu_vtx_z",&_nu_vtx_z,"nu_vtx_z/D");
  _tree->Branch("_x_cross",&_x_cross,"x_cross/D");
  _tree->Branch("_y_cross",&_y_cross,"y_cross/D");
  _tree->Branch("_z_cross",&_z_cross,"z_cross/D");
  _tree->Branch("_cross_pdg",&_cross_pdg,"_cross_pdg/I");
  _tree->Branch("_cross_momentum",&_cross_momentum,"_cross_momentum/D");
  _tree->Branch("_cross_xstart",&_cross_xstart,"_cross_xstart/D");
  _tree->Branch("_cross_ystart",&_cross_ystart,"_cross_ystart/D");
  _tree->Branch("_cross_zstart",&_cross_zstart,"_cross_zstart/D");
  _tree->Branch("_cross_xend",&_cross_xend,"_cross_xend/D");
  _tree->Branch("_cross_yend",&_cross_yend,"_cross_yend/D");
  _tree->Branch("_cross_zend",&_cross_zend,"_cross_zend/D");

}

void G4VetoSignalImpact::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(G4VetoSignalImpact)
