#ifndef MICHELCLUSTER_CXX
#define MICHELCLUSTER_CXX

#include "MichelConstants.h"
#include "MichelCluster.h"
#include "MichelException.h"
#include <cmath>
namespace michel {

  MichelCluster::MichelCluster(size_t min_nhits, double d_cutoff)
    : _min_nhits ( min_nhits )
    , _d_cutoff  ( d_cutoff  )
  { 
    _verbosity = msg::kNORMAL;
}

  MichelCluster::MichelCluster(const std::vector<HitPt>& hits,
			       size_t min_nhits,
			       double d_cutoff)
    : _hits(hits)
    , _min_nhits ( min_nhits )
    , _d_cutoff  ( d_cutoff  )
  {
    _verbosity = msg::kNORMAL;
    ProcessHits();
  }

  MichelCluster::MichelCluster(std::vector<HitPt>&& hits,
			       size_t min_nhits,
			       double d_cutoff)
    : _hits(std::move(hits))
    , _min_nhits ( min_nhits )
    , _d_cutoff  ( d_cutoff  )
  {
    _verbosity = msg::kNORMAL;
    ProcessHits();
  }
  /*
  MichelCluster::MichelCluster(const MichelCluster&& rhs)
    : _hits      ( std::move( rhs._hits     ) )
    , _start     ( std::move( rhs._start    ) )
    , _end       ( std::move( rhs._end      ) )
    , _ordered_pts ( std::move( rhs._ordered_pts ) )
    , _ds_v      ( std::move( rhs._ds_v     ) )
    , _s_v       ( std::move( rhs._s_v      ) )
    , _boundary  ( rhs._boundary              )
    , _michel    ( std::move( rhs._michel   ) )
    , _chi2_v    ( std::move( rhs._chi2_v   ) )
    , _t_mean_v  ( std::move( rhs._t_mean_v ) )
    , _t_dqds_v  ( std::move( rhs._t_dqds_v ) )
    , _verbosity ( rhs._verbosity )
    , _min_nhits ( rhs._min_nhits )
    , _d_cutoff  ( rhs._d_cutoff  )
  {}
  */
  void MichelCluster::Dump() const
  {
    std::cout << "\t\n==start dump==\n";
    std::cout << "A cluster with " << _hits.size() << " hits \n";
    std::cout << "The start point is at (" << _start._w << "," << _start._t << ")\n";
    std::cout << "The end point is at ("   << _end._w << ","   << _end._t   << ")\n";
    std::cout << _ordered_pts.size() << " of the hits are ordered and nearby ";
    //std::cout << "I may or may not have a michel, do i? " << _has_michel << "\n";
    /*
    if(_has_michel) 
      std::cout << "ok I do, it is at a distance : " << _michel_dist << " from one of the hits in my cham\n";
    */
    std::cout << "\t\n==end dump==\n";
  }

  void MichelCluster::SetHits(const std::vector<HitPt>& hits)
  {
    _hits = hits;
    ProcessHits();
  }

  void MichelCluster::SetHits(std::vector<HitPt>&& hits)
  {
    std::swap(_hits,hits);
    ProcessHits();
  }

  void MichelCluster::ProcessHits()
  {
    // Check if hit count
    if(_hits.size() < _min_nhits) {
      std::cerr << "\033[93m[ERROR]\033[00m "
		<< "<<" << __FUNCTION__ << ">> " 
		<< "Number of hit ("<< _hits.size() 
		<< ") is smaller than required (" 
		<< _min_nhits << ")" << std::endl;
      throw MichelException();
    }

    // Find min/max hits in terms of wire position
    size_t min_index = kINVALID_SIZE;
    size_t max_index = kINVALID_SIZE;
    double min_wire  = kMAX_DOUBLE;
    double max_wire  = kMIN_DOUBLE;
    for(size_t i=0; i<_hits.size(); ++i) {
      if(_hits[i]._w < min_wire) {
	min_wire  = _hits[i]._w;
	min_index = i;
      }
      if(_hits[i]._w > max_wire) {
	max_wire  = _hits[i]._w;
	max_index = i;
      }
    }
    // Check index validity
    if(min_index == kINVALID_SIZE || max_index == kINVALID_SIZE) {
      std::cerr << "\033[93m[ERROR]\033[00m <<" << __FUNCTION__ << ">> "
		<< "Did not find valid edge points..."
		<< std::endl;
      throw MichelException();
    }
    //order the points right => left
    OrderPoints(min_index, _ordered_pts, _ds_v, _s_v);

    //order the points left => right
    std::vector<size_t> lr_ordered_pts;
    std::vector<double> lr_ds_v, lr_s_v;
    OrderPoints(max_index, lr_ordered_pts, lr_ds_v, lr_s_v);

    // Make deciison: we take a longer cluster
    if(_ordered_pts.size() < lr_ordered_pts.size()) {
      _ordered_pts = lr_ordered_pts;
      _ds_v        = lr_ds_v;
      _s_v         = lr_s_v;
    }
    _start = _hits[ _ordered_pts.front() ];
    _end   = _hits[ _ordered_pts.back()  ];

    std::vector<HitPt> ordered_hits;
    ordered_hits.reserve(_ordered_pts.size());
    for(auto const& hit_index : _ordered_pts)
      ordered_hits.emplace_back(_hits[hit_index]);
    std::swap(ordered_hits,_hits);
    
    for(size_t i=0; i<_ordered_pts.size(); ++i)
      _ordered_pts[i]=i;
    
  }

  void MichelCluster::OrderPoints(size_t start_index,
				  std::vector<size_t>& ordered_index_v,
				  std::vector<double>& ds_v,
				  std::vector<double>& s_v)
  {
    // Clear result holder
    ordered_index_v.clear();
    s_v.clear();
    
    if(_hits.empty()) return;

    if(start_index >= _hits.size()) {
      std::cerr << "\033[93m[ERROR]\033[00m <<" << __FUNCTION__ << ">> "
		<< "Start index (" << start_index << ") is >= hit length "
		<< "(" << _hits.size() << ")" << std::endl;
      throw MichelException();
    }

    // Result holder, reserve max possible entries
    ordered_index_v.reserve(_hits.size());
    ordered_index_v.push_back(start_index);
    s_v.reserve  ( _hits.size() );
    ds_v.reserve ( _hits.size() );

    // Distance vector has the same length as the points vector
    s_v.push_back(0);
    
    // Marker vector to mark already-used-index
    std::vector<bool> used_v(_hits.size(),false);
    used_v[start_index]=true;
    while(ordered_index_v.size() < _hits.size()) {

      double min_dist  = kINVALID_DOUBLE;
      size_t min_index = kINVALID_SIZE;
      // Linear search: slow but robust ... think of a better robust way to replace
      for(size_t h_index=0; h_index<_hits.size(); ++h_index) {

	// If used skip
	if(used_v[h_index]) continue;

	// Take a reference of this hit & last hit
	auto const& this_step = _hits[h_index];
	auto const& last_step = _hits[ordered_index_v.back()];

	if( std::abs(this_step._w - last_step._w) > min_dist ) continue;
	if( std::abs(this_step._t - last_step._t) > min_dist ) continue;
	
	// Compute distance
	double sq_dist = pow(this_step._w - last_step._w,2) + pow(this_step._t - last_step._t,2);

	// If bigger than min distance, ignore this
	if(sq_dist > min_dist) continue;
	// Else register as the local min point
	min_dist  = sq_dist;
	min_index = h_index;
      }
      // If min_dist is above the cut-off, break
      if(min_dist > _d_cutoff) break;

      // Else this is a good index. Register and move on
      ordered_index_v.push_back(min_index);
      ds_v.push_back ( sqrt(min_dist)           );
      s_v.push_back  ( s_v.back() + ds_v.back() );
      used_v[min_index] = true;
    }
    // Verbosity report
    if( _verbosity <= msg::kINFO ) {
      // INFO level prints out where it starts & # points
      std::cout << "Ordered from index " << start_index
		<< " ... found " << ordered_index_v.size()
		<< " points!" <<std::endl;
      // DEBUG level prints out all indeces ... we never understand this anyway
      if(_verbosity <= msg::kDEBUG ) {
	for(size_t i=0; i<ordered_index_v.size(); ++i) {
	  std::cout << ordered_index_v[i] << " ";
	  if(i%8==0) std::cout<<std::endl;
	}
	std::cout<<std::endl;
      }
    }
  }

  MichelCluster MichelCluster::operator+(const MichelCluster& rhs) const
  {
    auto const& hits_lhs = (*this)._hits;
    auto const& hits_rhs = rhs._hits;

    std::vector<HitPt> hits;
    hits.reserve(hits_lhs.size() + hits_rhs.size());

    for(auto const& h : hits_lhs) hits.push_back(h);
    for(auto const& h : hits_rhs) hits.push_back(h);

    MichelCluster res(std::move(hits), _min_nhits, _d_cutoff);
    res.ProcessHits();
    return res;
  }

  MichelCluster& MichelCluster::operator+=(const MichelCluster& rhs)
  {
    _hits.reserve(_hits.size() + rhs._hits.size());

    for(auto const& h : rhs._hits) _hits.push_back(h);

    ProcessHits();

    return (*this);
  }

  const HitPt& MichelCluster::ClosestHit(const HitPt& ref)
  {
    double min_dist  = kMAX_DOUBLE;
    size_t min_index = kINVALID_SIZE;
    for(size_t i=0; i<_hits.size(); ++i) {

      auto const& h = _hits[i];

      double dist = h.SqDist(ref);

      if(dist < min_dist) {
	min_dist = dist;
	min_index = i;
      }
    }

    if(min_index == kINVALID_SIZE) {
      std::cerr << "\033[93m[ERROR]\033[00m <<" << __FUNCTION__ << ">> "
		<< "Did not find the closest hit point to the reference (ref ID="
		<< ref._id << ")" << std::endl;
      throw MichelException();
    }

    return _hits[min_index];
  }
  
  /*
int ClusterYPlane::match(const TVector2& michel_loc) {
  
  ///in cosmics
  if(_has_michel)
    return 2; //2 I already have a michel...
  
  auto dist = 999.0;
  auto idx  = size_t{0};
  
  for(size_t i = 0; i < _ordered_pts.size(); ++i) {
    auto d = distance(_ahits[_ordered_pts[i]].vec,michel_loc);
    if( d < dist) { dist = d; idx = i; }
  } 
  
  //distance to this michel is dist... , what is the cut off for this??
  if(dist < 4 ) { //? is this reasonable...
    _michel_location = idx;
    _michel_dist     = dist;
    _has_michel      = true;
  }
  else { return 0; } //michel isn't right for me :(
  
  return 1; //1 is matched with accepted parameters
}
  */
}
#endif

