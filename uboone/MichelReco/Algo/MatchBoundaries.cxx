#ifndef MICHELCLUSTER_MATCHQBOUNDARIES_CXX
#define MICHELCLUSTER_MATCHQBOUNDARIES_CXX

#include "MatchBoundaries.h"
#include "uboone/MichelReco/Fmwk/MichelException.h"
#include "uboone/MichelReco/Fmwk/ClusterVectorCalculator.h"
#include <cmath>
#include <cstdlib>
#include <sstream>
namespace michel {
  
  void MatchBoundaries::EventReset()
  {}
  
  bool MatchBoundaries::ProcessCluster(MichelCluster& cluster,
				       const std::vector<HitPt>& hits)
  { 

    if (hits.size() == 0) return false;
    
    /// call instance of "ClusterVectorCalculator"
    /// this class has a bunch of utility functions
    /// to calculate stuff based on the vecor of
    /// hits for a cluster
    ClusterVectorCalculator _clusterCalc;
    
    std::vector<double> truncated_mean;
    std::vector<double> truncated_dqds;
    std::vector<double> covariance; 
    std::vector<double> slope; 
    
    truncated_mean.reserve(cluster._ordered_pts.size());
    truncated_dqds.reserve(cluster._ordered_pts.size());
    covariance.reserve    (cluster._ordered_pts.size());
    slope.reserve         (cluster._ordered_pts.size());
    
    // //hardcoded for now will become configurable
    // double _n_window_size = 15;
    // double _p_above       = 0.25;
    // int    _window_cutoff = 3;

    //do truncated mean
    truncated_mean = _clusterCalc.calc_smooth_mean(cluster,
						   _n_window_size,
						   _window_cutoff,
						   _p_above);
    
    if(truncated_mean.size() < _edgefix) {
      std::stringstream ss;
      ss << std::endl
	 << "\tUnable to fix edges on truncated mean, edgefix size: " << _edgefix
	 << "\t and truncated_mean.size(): " << truncated_mean.size();
      Print(msg::kERROR,__FUNCTION__,ss.str());
      return false;
    }
    
    for(size_t i = 0 ; i < _edgefix; ++i) {
      truncated_mean.at(i) = truncated_mean[_edgefix];
      truncated_mean.at(truncated_mean.size() - i - 1) = truncated_mean[truncated_mean.size() - _edgefix];
    }
    
    //Directionality considerations
    int dir_window = _covariance_window;
    covariance     = _clusterCalc.calc_covariance(cluster._hits,dir_window);
    slope          = _clusterCalc.calc_slope     (cluster._hits,dir_window);
    
    // must be odd, currently has no setter,
    // sorry that this method has no info on it, ask vic
    int s = 3; 
    truncated_dqds = _clusterCalc.calc_smooth_derive(cluster._s_v,truncated_mean,s);
        
    //Lets play with truncated mean shaving...
    if(_verbosity <= msg::kINFO) {
      std::stringstream ss;
      ss << std::endl
	 << "\t\tIn MatchBoundaries" << std::endl
	 << "\tI have " << truncated_mean.size() << " truncated mean size" << std::endl
	 << "\twith   " << truncated_dqds.size() << " derivative points." << std::endl
	 << "\tMy incoming cluster has " << cluster._hits.size() << " hits in it...";
      Print(msg::kINFO,__FUNCTION__,ss.str());
    }
    
    //With this new information, calculate the boundary point between possible muon end and michel start
    
    size_t candidate_loc     = _clusterCalc.find_max(truncated_mean);
    size_t dqdscandidate_loc = _clusterCalc.find_min(truncated_dqds); 

    std::swap(cluster._t_mean_v,truncated_mean);
    std::swap(cluster._t_dqds_v,truncated_dqds);
    
    if((candidate_loc     >= cluster._hits.size()))
      return false;
    
    if((dqdscandidate_loc >= cluster._hits.size()))
      return false;
    
    if(std::abs(int(dqdscandidate_loc) - int(candidate_loc)) > _maxDistance)
      return false;
    

    size_t right = cluster._ordered_pts.size() - 1 - candidate_loc;
    size_t left  = candidate_loc;
    
    size_t iMin = 0;
    size_t iMax = 0;
    
    if(right >= _maxDistance) iMax  = _maxDistance   + candidate_loc;
    if(left  >= _maxDistance) iMin  = candidate_loc - _maxDistance;

    if(right < _maxDistance)  iMax  = cluster._hits.size() - 1;
    if(left  < _maxDistance)  iMin  = 0;

    // holder for hit with largest charge -> this will identify the location
    // of the hit that defines the michel boundary
    auto k   = 0.0;
    size_t idx = 0;
    
    for(size_t w = iMin; w <= iMax; ++w) {
      auto c = cluster._hits[cluster._ordered_pts[w]]._q;
      // if this hit has more charge than any other
      if(c > k) { k = c; idx = w; }
    }

    // move on only if avg. covariance for hits near the "start" is less than some value
    // idx is the position of the reconstructed michel start
    double avg_covariance = 0;
    int counts = 0;

    size_t start = 0;
    if (idx >= 3) start = idx-3;
    
    for (size_t i = start; i < (idx+4); i++){
      // make sure we fall in the vector's range
      if ( i < covariance.size() ){
	avg_covariance += std::abs(covariance[i]);
	counts += 1;
      }
    }
    // make sure we have at least 1 point!
    if (avg_covariance == 0)
      return false;

    // if so, check that the average covariance is below
    // the max allowed covariance
    avg_covariance /= counts;
    if (avg_covariance > _maxCovarianceAtStart)
      return false;
    
    std::swap(cluster._chi2_v,covariance);
    std::swap(cluster._dirs_v,slope);

    cluster._boundary = cluster._ordered_pts[idx];
    return true;
  }
  
}
#endif
