#include "ubcrt/CRTXSEC/CRTAnaFun.hh"
#include "cetlib_except/exception.h"


int crtana::auxfunc::test(int a){
  
  int b = a + 5731;//tururu
  return b;
}
/*
int crtana::auxfunc::test2(int a){
  
  int b = a + 5731;//tururu
  return b;
}*/
/*std::vector<crt::CRTHit> crtana::auxfunc::getCRTHitVector(art::Event const & evt, std::string data_label_hits_){
  art::Handle< std::vector<crt::CRTHit> > rawHandle_hits;
  evt.getByLabel(data_label_hits_, rawHandle_hits);
  if(!rawHandle_hits.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has " << "\033[31m" << "zero" << "\033[0m"
              << " CRTHits " << " in module " << data_label_hits_ << std::endl;
    std::cout << std::endl;
  }
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle_hits);
  if(rawHandle_hits.isValid()){
    //get better access to the data               
    std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle_hits);
    if(rawHandle_hits.isValid()){
      std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
                << ", event " << evt.event() << " has " << "\033[32m" << CRTHitCollection.size() << "\033[0m"
                << " CRTHits " << " in module " << data_label_hits_ << std::endl;
     //return;
    }
    
  }
    return CRTHitCollection;
  
  
}*/

  void crtana::auxfunc::GetPFParticleIdMap(const crtana::PFParticleHandle &pfParticleHandle, crtana::PFParticleIdMap &pfParticleMap)
  {
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
      {
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second)
          {
              throw cet::exception("crt_ana") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
          }
      }
  }


  void crtana::auxfunc::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles)
  {
    std::cout << "Number of pfParticles: " << pfParticleMap.size() << std::endl;
    int neutrino_candidate_counter = 0;
    //int cosmic_counter = 0;
    //int nu_counter = 0;
      for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)
      {
          const art::Ptr<recob::PFParticle> pParticle(it->second);

          // Only look for primary particles
          if (!pParticle->IsPrimary()) continue;

          // Check if this particle is identified as the neutrino
          const int pdg(pParticle->PdgCode());
          const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);

          // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
          if (!isNeutrino)
          {
              crParticles.push_back(pParticle);
              continue;
          }
          neutrino_candidate_counter++;

          // ATTN. We are filling nuParticles under the assumption that there is only one reconstructed neutrino identified per event.
          //       If this is not the case please handle accordingly
          if (!nuParticles.empty())
          {
              //throw cet::exception("crt_ana") << "  This event contains multiple reconstructed neutrinos!";
          }
        std::cout << "neutrino candidates: " << neutrino_candidate_counter << " has : ";
          int daughter_counter = 0;
          // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
          for (const size_t daughterId : pParticle->Daughters())
          {
              if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                  throw cet::exception("crt_ana") << "  Invalid PFParticle collection!";

              nuParticles.push_back(pfParticleMap.at(daughterId));
            daughter_counter++;
          }
        std::cout << daughter_counter << " Daughters" << std::endl;
      }
    std::cout << "Number pf neutrino candidates daugthers: " << nuParticles.size() << std::endl;
    std::cout << "Number pf primary cosmic candidates: " << crParticles.size() << std::endl;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void crtana::auxfunc::CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers, std::string m_trackLabel, std::string m_showerLabel)
  {
      // Get the associations between PFParticles and tracks/showers from the event

      //std::string m_trackLabel = "";           ///< The label for the track producer from PFParticles
      //std::string m_showerLabel = "";

      art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, evt, m_trackLabel);
      art::FindManyP< recob::Shower > pfPartToShowerAssoc(pfParticleHandle, evt, m_showerLabel);

      for (const art::Ptr<recob::PFParticle> &pParticle : particles)
      {
          const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
          const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));
          const unsigned int nTracks(associatedTracks.size());
          const unsigned int nShowers(associatedShowers.size());

          // Check if the PFParticle has no associated tracks or showers
          if (nTracks == 0 && nShowers == 0)
          {
              mf::LogDebug("crt_ana") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << std::endl;
              continue;
          }

          // Check if there is an associated track
          if (nTracks == 1 && nShowers == 0)
          {
              tracks.push_back(associatedTracks.front());
              continue;
          }

          // Check if there is an associated shower
          if (nTracks == 0 && nShowers == 1)
          {
              showers.push_back(associatedShowers.front());
              continue;
          }

          throw cet::exception("crt_ana") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();
      }
  }

double crtana::auxfunc::TpcTrack_match_CrtHit(const art::Ptr<recob::Track> &ptracks, std::vector<crt::CRTHit> const& CRTHitCollection, std::vector<double> &dist_return, double crthitmatch_){
  std::vector<TVector3> sorted_trk;
  recob::Track my_track = *ptracks.get();
  SortTrackPoints(my_track,sorted_trk);
  std::vector<double> dist;

  auto const& top    = sorted_trk.at(0);
  auto & bottom = sorted_trk.at(sorted_trk.size() - 1);
  int counter = 1;
  while(bottom.X()==-999 && (sorted_trk.size()-counter) >=0) {
    bottom = sorted_trk.at(sorted_trk.size() - counter);
    counter++;
  }
  
  if(top.X()==-999 || top.Y()==-999 || top.Z()==-999 || bottom.X()==-999 || bottom.Y()==-999 || bottom.Z()==-999){
    std::cout<<"Track: "<<my_track.ID() << " has invalid track variables. Length: " << my_track.Length() << std::endl;
    //for(std::vector<int>::size_type i = 0; i != sorted_trk.size(); i++) {
    //  std::cout<<"TPC Track point: ("<<sorted_trk.at(i).X()<<","<<sorted_trk.at(i).Y()<<","<<sorted_trk.at(i).Z()<<")"<<std::endl;
    //}
    return 0;
  }

  double TT[] = {top.X(),top.Y(),top.Z()};
  double BT[] = {bottom.X(),bottom.Y(),bottom.Z()};
  //double VT[] = {BT[0] - TT[0], BT[1] - TT[1], BT[2] - TT[2]};
  double VT[] = {TT[0] - BT[0], TT[1] - BT[1], TT[2] - BT[2]};

  double TPCTheta = crt::auxfunctions::CalTheta(VT[0],VT[1],VT[2]);
  double TPCPhi = crt::auxfunctions::CalPhi(VT[0],VT[1],VT[2]);
  
  if(0){
    std::cout<<"TPC Track Director Vector: ("<<VT[0]<<","<<VT[1]<<","<<VT[2]<<")"<<std::endl;
    std::cout<<"TPC Theta: "<<TPCTheta<<std::endl;
    std::cout<<"TPC Phi: "<<TPCPhi<<std::endl;
  }
  int matching = 0;
  double d = 999;
  for(std::vector<int>::size_type i = 0; i != CRTHitCollection.size(); i++) {//A 
    crt::CRTHit my_CRTHit = CRTHitCollection[i];
    if(my_CRTHit.plane == 0 || my_CRTHit.plane == 3){ // Bottom or Top
      double p = ( my_CRTHit.y_pos-top.Y() ) / VT[1];
      double Z_prop = top.Z() + VT[2]*p;
      d = abs(my_CRTHit.z_pos - Z_prop);
      //std::cout<<"difference in bottom or top: "<<d<<std::endl;
      dist.push_back(d);
      //std::cout<<"difference in bottom or top: "<<d<<std::endl;
      if(d<crthitmatch_) matching = 1;
      //std::cout<<"difference in bottom or top: "<<d<<std::endl;
    }
    else if(my_CRTHit.plane == 1 || my_CRTHit.plane == 2){ //Feed-through
      double a = VT[1]/VT[2];
      double b = 1;
      double c = TT[1] - (TT[2]*VT[1] / VT[2] );
      d = abs(a * my_CRTHit.z_pos + b * my_CRTHit.y_pos + c) / sqrt(a*a + b*b);
      //std::cout<<"difference in feed or pipe: "<<d<<std::endl;
      dist.push_back(d);
      if(d<crthitmatch_) matching = 1;
    }
    /*if(d<45 && my_track.Length()>200){
      std::cout<<"TPC Track Director Vector: ("<<VT[0]<<","<<VT[1]<<","<<VT[2]<<")"<<std::endl;
      std::cout<<"TPC Track Beg Vector: ("<<TT[0]<<","<<TT[1]<<","<<TT[2]<<")"<<std::endl;
      std::cout<<"TPC Track End Vector: ("<<BT[0]<<","<<BT[1]<<","<<BT[2]<<")"<<std::endl;
      std::cout<<"TPC Theta: "<<TPCTheta<<std::endl;
      std::cout<<"TPC Phi: "<<TPCPhi<<std::endl;
      std::cout<<"CRT Hit Pos: ("<<my_CRTHit.x_pos<<","<<my_CRTHit.y_pos<<","<<my_CRTHit.z_pos<<")"<<std::endl;
  }*/
    
    
  }
  
  
  dist_return =  dist;
  
  return matching;
}

void crtana::auxfunc::SortTrackPoints(const recob::Track& track, std::vector<TVector3>& sorted_trk)
{

  sorted_trk.clear();

  auto const&N = track.NumberTrajectoryPoints();
  auto const&start = track.LocationAtPoint(0);
  auto const&end   = track.LocationAtPoint( N - 1 );

  if (start.Y() > end.Y()){
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( track.LocationAtPoint(i) );
  }

  else {
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( track.LocationAtPoint( N - i - 1) );
  }
}

