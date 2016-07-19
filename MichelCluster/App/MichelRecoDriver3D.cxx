#ifndef LARLITE_MICHELRECODRIVER3D_CXX
#define LARLITE_MICHELRECODRIVER3D_CXX

#include "MichelRecoDriver3D.h"

#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/simch.h"
#include "DataFormat/track.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"


namespace larlite {

  MichelRecoDriver3D::MichelRecoDriver3D()
    : _hit_tree(nullptr)
  {
    _name="MichelRecoDriver3D";
    _fout=0;
    _save_clusters=false;
    _Efield=0.5;
    // FIX ME: currently set only plane 2 reconstruction
    SetPlane(2);
  }

  bool MichelRecoDriver3D::initialize() {

    if(_producer.empty()) {
      print(msg::kERROR,__FUNCTION__,"Input data product producer name not specified!");
      throw std::exception();
    }

    _events_info.clear();
    
    _hit_tree = new TTree("_hit_tree","Hit Tree [all hits in event]");
    _hit_tree->Branch("_run"   , &_run    , "run/I");
    _hit_tree->Branch("_subrun", &_subrun , "subrun/I");
    _hit_tree->Branch("_event" , &_event  , "event/I");
    _hit_tree->Branch("_q_v" , "std::vector<double>" , &_q_v);
    _hit_tree->Branch("_w_v" , "std::vector<double>" , &_w_v);
    _hit_tree->Branch("_t_v" , "std::vector<double>" , &_t_v);
    _hit_tree->Branch("_p_v" , "std::vector<double>" , &_p_v);

    _mgr.Initialize();

    _event_time = 0;
    _event_ctr  = 0;

    return true;
  }
  
  bool MichelRecoDriver3D::analyze(storage_manager* storage) {

    // load GeometryHelper utility
    auto geomHelper = ::larutil::GeometryHelper::GetME();

    _event_watch.Start();

    // get ID information
    _run    = storage->get_data<event_cluster>(_producer)->run();
    _subrun = storage->get_data<event_cluster>(_producer)->subrun();
    _event  = storage->get_data<event_cluster>(_producer)->event_id();

    // grab reconstructed tracks
    auto ev_track = storage->get_data<event_track>(_producer);

    // create output tracks and clusters
    auto out_track = storage->get_data<event_track>("output");
    auto out_clus  = storage->get_data<event_cluster>("michel");
    // create association object for clusters
    auto out_clus_hit_ass_v_v = storage->get_data<event_ass>(out_clus->name());

    std::vector<std::vector<unsigned int> > michel_clus_hit_ass_v_v; // michel cluster -> hit

    if(!ev_track || ev_track->empty()) {
      std::cout<<"No tracks found. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }
    
    // grab the hits associated with the reconstructed tracks
    larlite::event_hit *ev_hit = nullptr;
    auto const& trk_hit_ass_v_v = storage->find_one_ass( ev_track->id(), ev_hit, ev_track->name() );

    // If ev_hit is null, failed to find data or assocaition (shouldn't happen)
    if(!ev_hit || ev_hit->empty()) {
      std::cout << "No hit found associated to 3D tracks. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }


    storage->set_id(_run,_subrun,_event);

    // 1st thing to do is loop through tracks and create GeoAlgo track objects.
    _geoTrj_v.clear();
    for (size_t i=0; i < ev_track->size(); i++){
      // create geotrack from track
      auto const& trk = ev_track->at(i);
      ::geoalgo::Trajectory_t trj;
      if (trk.NumberTrajectoryPoints() < 2) continue;
      for (size_t n=0; n < trk.NumberTrajectoryPoints(); n++){
	auto const& pt = trk.LocationAtPoint(n);
	trj.push_back ( pt );
      }// for all traj points
      _geoTrj_v.push_back(trj);
    }// for all tracks

    std::cout << " added " << _geoTrj_v.size() << " tracks" << std::endl;


    
    // 1st step: check for 3D kink topology in tracks
    std::vector< std::pair<int,int> > trj_pairs_v;
    std::vector<geoalgo::Point_t>   PoCA_v;
    SearchTrajectoryPairs(trj_pairs_v,PoCA_v);
    std::cout << "found " << trj_pairs_v.size() << " good matches" << std::endl;
    for (size_t t=0; t < trj_pairs_v.size(); t++){
      auto const& trj_pair = trj_pairs_v[t];
      out_track->emplace_back( ev_track->at( trj_pair.first  ) );
      out_track->emplace_back( ev_track->at( trj_pair.second ) );

      // get longest track of the two
      double len1 = ev_track->at( trj_pair.first  ).Length();
      double len2 = ev_track->at( trj_pair.second ).Length();
      
      std::vector<unsigned int> trk_hit_indices;
      if (len1 < len2){
	for (auto const& hit_idx : trk_hit_ass_v_v[ trj_pair.first ] ){
	  auto const& hit = ev_hit->at(hit_idx);
	  if (hit.WireID().Plane == 2)
	    trk_hit_indices.push_back( hit_idx );
	}
      }
      else{
	for (auto const& hit_idx : trk_hit_ass_v_v[ trj_pair.second ] ){
	  auto const& hit = ev_hit->at(hit_idx);
	  if (hit.WireID().Plane == 2)
	    trk_hit_indices.push_back( hit_idx );
	}
      }
      
      std::cout << "Shortest track has " << trk_hit_indices.size() << " hits " << std::endl;
      if (trk_hit_indices.size() == 0)
	continue;

      michel_clus_hit_ass_v_v.push_back( trk_hit_indices );

      // create output larlite cluster
      double *xyz = new double[3];
      xyz[0] = PoCA_v[t][0];
      xyz[1] = PoCA_v[t][1];
      xyz[2] = PoCA_v[t][2];
      auto const& projection2D = geomHelper->Point_3Dto2D(xyz, 2);
      int wire = (int) ev_hit->at( trk_hit_indices[0] ).WireID().Wire;
      int time = (int) ev_hit->at( trk_hit_indices[0] ).PeakTime();
      larlite::cluster clus;
      clus.set_start_wire( wire, 0. );
      clus.set_start_tick( time, 0. );
      out_clus->emplace_back( clus );
      // get shortest track of the two. This is the Michel
      
      
    }

    out_clus_hit_ass_v_v->set_association(out_clus->id(),product_id(data::kHit,ev_hit->name()),michel_clus_hit_ass_v_v);

    return true;
  }


  bool MichelRecoDriver3D::finalize() {

    std::cout << "time/event = " << _event_time/_event_ctr * 1.e6 << std::endl;


    auto ts = ::larutil::TimeService::GetME();
    std::cout << "TPCTick -> TDC: 1000 -> " << ts->TPCTick2TDC(1000) << std::endl;
    std::cout << "TPCTick -> TDC: 2000 -> " << ts->TPCTick2TDC(2000) << std::endl;

    _fout->cd();
    _mgr.Finalize(_fout);
    if (_hit_tree)
      _hit_tree->Write();
    return true;
  }


  void MichelRecoDriver3D::SearchTrajectoryPairs(std::vector< std::pair<int,int> >& trj_pairs_v,
						 std::vector<geoalgo::Point_t>& PoCA_v) {
    
    trj_pairs_v.clear();
    PoCA_v.clear();

    for (size_t n=0; n < _geoTrj_v.size(); n++){

      auto const& trj1 = _geoTrj_v[n];

      if (trj1.Length() < 4) continue;

      for (size_t m=n+1; m < _geoTrj_v.size(); m++){

	auto const& trj2 = _geoTrj_v[m];

	if (trj2.Length() < 4) continue;

	geoalgo::Point_t pt1;
	geoalgo::Point_t pt2;

	double sqDistMin = 1000.;
	double sqDist1 = trj1.at(0).Dist( trj2.at(0) );
	if (sqDist1 < sqDistMin) { sqDistMin = sqDist1; pt1 = trj1.at(0); pt2 = trj2.at(0); }
	double sqDist2 = trj1.at( trj1.size()-1 ).Dist( trj2.at(0) );
	if (sqDist2 < sqDistMin) { sqDistMin = sqDist2; pt1 = trj1.at( trj1.size()-1 ); pt2 = trj2.at(0); }
	double sqDist3 = trj1.at(0).Dist( trj2.at( trj2.size()-1 ) );
	if (sqDist1 < sqDistMin) { sqDistMin = sqDist3; pt1 = trj1.at(0); pt2 = trj2.at( trj2.size()-1 ); }
	double sqDist4 = trj1.at( trj1.size()-1 ).Dist( trj2.at( trj2.size()-1 ) );
	if (sqDist1 < sqDistMin) { sqDistMin = sqDist4; pt1 = trj1.at( trj1.size()-1 ); pt2 = trj2.at( trj2.size()-1 ); }
	
	if (sqDistMin < 5){
	  PoCA_v.push_back( ( pt1 + pt2 ) / 2 );
	  trj_pairs_v.push_back( std::make_pair(n,m) );
	}

      }// for all trajectories, loop 2
    }// for all trajectories, loop 1

    return;

  }

}
#endif
