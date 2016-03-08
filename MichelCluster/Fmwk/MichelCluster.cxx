#ifndef MICHELCLUSTER_CXX
#define MICHELCLUSTER_CXX

#include "MichelConstants.h"
#include "MichelCluster.h"
#include "MichelException.h"
#include "ClusterVectorCalculator.h"
#include "UtilFunc.h"
#include <cmath>
#include <sstream>
namespace michel {

  std::string MichelCluster::Diff(const MichelCluster& rhs) const
  {
    std::string msg;
    if(_start != rhs._start)
      msg += "  HitPt _start changed\n";

    if(_end != rhs._end)
      msg += "  HitPt _end changed\n";

    if(_forward != rhs._forward)
      msg += "  bool _forward changed\n";

    if(_boundary != rhs._boundary)
      msg += "  HitIdx_t _boundary changed\n";

    if(_min_nhits != rhs._min_nhits)
      msg += "  size_t _min_nhits changed\n";

    if(_d_cutoff != rhs._d_cutoff)
      msg += "  double _d_cutoff changed\n";

    if(_has_michel != rhs._has_michel)
      msg += "  bool _has_michel changed\n";

    if(!IsSame(_hits,rhs._hits))
      msg += "    std::vector<HitPt> _hits changed\n";

    if(!IsSame(_ds_v,rhs._ds_v))
      msg += "    std::vector<double> _ds_v changed\n";

    if(!IsSame(_s_v,rhs._s_v))
      msg += "    std::vector<double> _s_v changed\n";

    if(!IsSame(_chi2_v,rhs._chi2_v))
      msg += "    std::vector<double> _chi2_v changed\n";
    
    if(!IsSame(_t_mean_v,rhs._t_mean_v))
      msg += "    std::vector<double> _t_mean_v changed\n";

    if(!IsSame(_t_dqds_v,rhs._t_dqds_v))
      msg += "    std::vector<double> _t_dqds_v changed\n";

    if(!IsSame(_dirs_v,rhs._dirs_v))
      msg += "    std::vector<double> _dirs_v changed\n";
    
    msg += _michel.Diff(rhs._michel);

    return msg;

  }

  MichelCluster::MichelCluster(size_t min_nhits, double d_cutoff)
    : _min_nhits ( min_nhits     )
    , _d_cutoff  ( d_cutoff      )
    , _id        ( kINVALID_SIZE ) 
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
    try{ ProcessHits(); }
    catch (MichelException exception) { _hits.clear(); }
  }

  MichelCluster::MichelCluster(std::vector<HitPt>&& hits,
			       size_t min_nhits,
			       double d_cutoff)
    : _hits(std::move(hits))
    , _min_nhits ( min_nhits )
    , _d_cutoff  ( d_cutoff  )
  {
    _verbosity = msg::kNORMAL;
    try { ProcessHits(); }
    catch (MichelException exception) { _hits.clear(); }
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
    std::stringstream ss;
    ss << "\t\n==start dump==" << std::endl
       << "A cluster with " << _hits.size() << " hits " << std::endl
       << "The start point is at (" << _start._w << "," << _start._t << ")" << std::endl
       << "The end point is at ("   << _end._w << ","   << _end._t   << ")" << std::endl
       << _ordered_pts.size() << " of the hits are ordered and nearby "
      //<< "I may or may not have a michel, do i? " << _has_michel << "" << std::endl;
      /*
	if(_has_michel) 
	<< "ok I do, it is at a distance : " << _michel_dist << " from one of the hits in my cham" << std::endl;
      */
       << "\t\n==end dump==" << std::endl;
    Print(msg::kNORMAL,__FUNCTION__,ss.str());
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
      std::stringstream ss; 
      ss << "<<" << __FUNCTION__ << ">> " 
	 << "Number of hit ("<< _hits.size() 
	 << ") is smaller than required (" 
	 << _min_nhits << ")";
      Print(msg::kEXCEPTION,__FUNCTION__,ss.str());
    }

    _all_hits = _hits;

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
    if(min_index == kINVALID_SIZE || max_index == kINVALID_SIZE){
      std::stringstream ss; 
      ss << "Did not find valid edge points. Hit list size is " << _hits.size();
      Print(msg::kEXCEPTION,__FUNCTION__,ss.str());
    }

    ClusterVectorCalculator _clusCalc;

    //order the points right => left
    _clusCalc.OrderPoints(min_index, _hits, _d_cutoff, _ordered_pts, _ds_v, _s_v);

    //order the points left => right
    std::vector<size_t> lr_ordered_pts;
    std::vector<double> lr_ds_v, lr_s_v;
    _clusCalc.OrderPoints(max_index, _hits, _d_cutoff, lr_ordered_pts, lr_ds_v, lr_s_v);

    // Make deciison: we take a longer cluster
    if(_ordered_pts.size() < lr_ordered_pts.size()) {
      _ordered_pts = lr_ordered_pts;
      _ds_v        = lr_ds_v;
      _s_v         = lr_s_v;
    }
    _start = _hits[ _ordered_pts.front() ];
    _end   = _hits[ _ordered_pts.back()  ];

    // now that we have a list of ordered hits
    // swap the current _hits list which
    // contains all hits in the cluster
    // with the list containing only the
    // ordered hit-list

    std::vector<HitPt> ordered_hits;
    ordered_hits.reserve(_ordered_pts.size());
    for(auto const& hit_index : _ordered_pts)
      ordered_hits.emplace_back(_hits[hit_index]);
    std::swap(ordered_hits,_hits);
    
    for(size_t i=0; i<_ordered_pts.size(); ++i)
      _ordered_pts[i]=i;
    
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

    res._id = this->_id;
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
      std::stringstream ss;
      ss << "Did not find the closest hit point to the reference (ref ID="
	 << ref._id << ")";
      Print(msg::kEXCEPTION,__FUNCTION__,ss.str());
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

