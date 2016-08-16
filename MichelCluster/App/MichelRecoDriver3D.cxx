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

    // grab clusters associated to tracks
    larlite::event_cluster *ev_clus = nullptr;
    auto const& trk_clus_ass_v_v = storage->find_one_ass( ev_track->id(), ev_clus, ev_track->name() );
    
    // grab the hits associated with the reconstructed tracks
    larlite::event_hit *ev_hit = nullptr;
    auto const& trk_hit_ass_v_v = storage->find_one_ass( ev_clus->id(), ev_hit, ev_clus->name() );
    
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

    // If ev_hit is null, failed to find data or assocaition (shouldn't happen)
    if(!ev_clus || ev_clus->empty()) {
      std::cout << "No hit found associated to 3D tracks. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }
    
    // If ev_hit is null, failed to find data or assocaition (shouldn't happen)
    if(!ev_hit || ev_hit->empty()) {
      std::cout << "No hit found associated to 3D tracks. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }

    storage->set_id(_run,_subrun,_event);

    // 1st thing to do is loop through tracks and create GeoAlgo track objects.
    _geoTrj_v.clear();
    _geoTrj_map.clear();
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
      _geoTrj_map[ _geoTrj_v.size() - 1 ] = i;
    }// for all tracks

    for (size_t i=0; i < _geoTrj_v.size(); i++){

      auto const& trj = _geoTrj_v[i];

      if (trj.Length() < 50)
	continue;

      ::geoalgo::Point_t muStop;
      if (trj.at(0)[1] > trj.at(trj.size()-1)[1])
	muStop = trj.at(trj.size()-1);
      else
	muStop = trj.at(0);

      if (muStop[1] < -90)
	continue;
      if ( (muStop[2] < 30) or (muStop[2] > 1000) )
	continue;
      if ( (muStop[0] < -20) or (muStop[0] > 290) )
	continue;

      // match trajectories to make sure there are no other ones nearby
      bool good = true;
      for (size_t j=0; j < _geoTrj_v.size(); j++){

	if (i == j) continue;

	auto const& trj2 = _geoTrj_v[j];

	if (trj2.Length() < 10)
	  continue;
	
	double dist = _geoAlgo.SqDist(muStop, trj2);
	if (dist < 30*30){
	  good = false;
	  break;
	}
	
      }// second loop of tracks
      
      if (good == false)
	continue;

      // take end point and convert to collection plane point
      double *xyz = new double[3];
      xyz[0] = muStop[0];
      xyz[1] = muStop[1];
      xyz[2] = muStop[2];
      auto const& projection2D = geomHelper->Point_3Dto2D(xyz, 2);

      int wire = (int) ( projection2D.w / geomHelper->WireToCm() );
      int time = (int) ( projection2D.t / geomHelper->TimeToCm() );

      larlite::cluster clus;
      clus.set_start_wire( wire, 0. );
      clus.set_start_tick( time, 0. );
      out_clus->emplace_back( clus );
      
      std::cout << "saving a michel @ [w,t] -> [" << wire << ", " << time << "]" << std::endl;
      
      std::vector<unsigned int> v;
      michel_clus_hit_ass_v_v.push_back( v );
      
    }// for all trajectories
    
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

}
#endif
