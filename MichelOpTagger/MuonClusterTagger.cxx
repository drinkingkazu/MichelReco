#ifndef LARLITE_MUONCLUSTERTAGGER_CXX
#define LARLITE_MUONCLUSTERTAGGER_CXX

#include "MuonClusterTagger.h"
#include "DataFormat/opflash.h"
#include "DataFormat/cluster.h"
#include "GeoAlgo/GeoAlgo.h"
#include "GeoAlgo/GeoLineSegment.h"
#include "DataFormat/hit.h"
#include "DataFormat/calorimetry.h"
#include "DataFormat/mctrack.h"
#include "LArUtil/GeometryUtilities.h"
#include "LArUtil/LArProperties.h"

namespace larlite {

  MuonClusterTagger::MuonClusterTagger()
    : _tree(nullptr)
    , _match_tree(nullptr)
  {
    _Efield = 0.5; // kV/cm
    _clusProducer = "muon";
    _use_mc = false;
    _use_y  = true;
    _name = "MuonClusterTagger";
    _fout = 0;
    return;
  }
  
  bool MuonClusterTagger::initialize() {

    // use instances of LArUtil and GeometryUtilities
    // for (w,t) -> (cm, cm) conversion
    // wire->cm
    _w2cm = larutil::GeometryUtilities::GetME()->WireToCm();
    // time->cm (accounting for different operating voltages)
    double driftVel = larutil::LArProperties::GetME()->DriftVelocity(_Efield,87); // [cm/us]
    // tick width in time
    double tickWidth = 0.5; // [us]
    _t2cm = tickWidth*driftVel;

    _tree = new TTree("flash_tree","");
    _tree->Branch("npe",&_npe,"npe/D");
    _tree->Branch("fx",&_flash_x,"fx/D");
    _tree->Branch("fy",&_flash_y,"fy/D");
    _tree->Branch("fz",&_flash_z,"fz/D");
    _tree->Branch("tx",&_tpc_x,"tx/D");
    _tree->Branch("ty",&_tpc_y,"ty/D");
    _tree->Branch("tz",&_tpc_z,"tz/D");
    _tree->Branch("ft",&_flash_time,"ft/D");
    _tree->Branch("mct",&_mc_time,"mct/D");
    _tree->Branch("mcx",&_mc_x,"mcx/D");
    _tree->Branch("mcy",&_mc_y,"mcy/D");
    _tree->Branch("mcz",&_mc_z,"mcz/D");
    _tree->Branch("mc_dx",&_mc_dx,"mc_dx/D");
    _tree->Branch("score",&_score,"score/D");

    _match_tree = new TTree("match_tree","");
    _match_tree->Branch("dz",&_dz,"dz/D");
    _match_tree->Branch("dy",&_dy,"dy/D");
    _match_tree->Branch("muon_t",&_muon_t,"muon_t/D");
    _match_tree->Branch("mimchel_t",&_michel_t,"michel_t/D");
    _match_tree->Branch("muon_pe",&_muon_pe,"muon_pe/D");
    _match_tree->Branch("mimchel_pe",&_michel_pe,"michel_pe/D");

    return true;
  }
  
  bool MuonClusterTagger::analyze(storage_manager* storage) {

    _mgr.Reset();

    // grab OpFlash data-product
    auto ev_flash = storage->get_data<event_opflash>("opflash");

    if(!ev_flash || ev_flash->empty()) {
      std::cout<<"No opflash found. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }

    // grab Cluster data-product
    auto ev_clus = storage->get_data<event_cluster>(_clusProducer);
    if(!ev_clus || ev_clus->empty()) {
      //std::cout<<"No cluster found. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }

    auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");

    // get hits associated with the clusters
    event_hit* ev_hit = nullptr;
    auto const& hit_ass_set = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());
    
    // If ev_hit is null, failed to find data or assocaition (shouldn't happen)
    if(!ev_hit || ev_hit->empty()){
      std::cout<<"No hits associated to cluster found. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }
    
    if (!_use_mc){
      // per cluster, make an optical track
      // Loop over clusters & add them to our algorithm manager
      for(auto const& hit_ass : hit_ass_set) {

	// QCluster object
	::flashana::QCluster_t tpc_obj;
	
	// Loop over hits and fill
	for (size_t i=0; i < hit_ass.size(); i++){
	  
	  auto const& h = (*ev_hit)[hit_ass[i]];  
	  
	  ::flashana::QPoint_t pt;
	  
	  pt.x = h.WireID().Wire * _w2cm;
	  pt.y = 0;
	  pt.z = h.PeakTime() * _t2cm;
	  pt.q = h.Integral();
	  
	  tpc_obj.emplace_back(pt);
	}// for all hits in the cluster
	// add the QCluster to the manager
	_mgr.Emplace(std::move(tpc_obj));
      }// for all clusters
      

    }// if we don't want to use MC
    // if we want to use MC information
    else{
      if (!ev_mctrack || ev_mctrack->empty()) return false;
      for(auto const& trk : *ev_mctrack) {
	
	::flashana::QCluster_t tpc_obj;
	
	if(trk.size()>=2) {
	  tpc_obj.reserve(trk.size()-1);
	  
	  for(size_t i=0; i < (trk.size()-1); ++i) {
	    
	    auto const& pt1 = trk[i].Position();
	    auto const& pt2 = trk[i+1].Position();
	    
	    ::flashana::QPoint_t pt;
	    
	    double dx = pt2[0] - pt1[0];
	    double dz = pt2[2] - pt1[2];
	    
	    pt.q = (trk[i].E() - trk[i+1].E());
	    pt.x = pt1[0] + dx/2.;
	    pt.y = 0;
	    pt.z = pt1[2] + dz/2.;

	    if (_use_y){
	      double dy = pt2[1] - pt1[1];
	      pt.y = pt1[1] + dy/2.;
	    }
	    
	  tpc_obj.emplace_back(pt);
	  }
	}
	_mgr.Emplace(std::move(tpc_obj));

      }
    }

    // now loop through optical flash objects
    for(auto const& flash : *ev_flash) {
      
      ::flashana::Flash_t f;
      f.x = f.x_err = 0;
      f.y = flash.YCenter();
      f.z = flash.ZCenter();
      f.y_err = flash.YWidth();
      f.z_err = flash.ZWidth();
      f.pe_v.reserve(32);
      for(unsigned int i=0; i<32; i++)
	f.pe_v.push_back(flash.PE(i));
      f.time = flash.Time();
      
      _mgr.Emplace(std::move(f));
    }// for all optical flashes
      
    auto const res = _mgr.Match();
    std::cout << "number of matches found: " << res.size() << std::endl;
    ::geoalgo::LineSegment line;
    ::geoalgo::Point_t pt(0,0,0);
    ::geoalgo::GeoAlgo geoalg;
    for(auto const& match : res) {
      auto const& flash = (*ev_flash)[match.flash_id];
      _flash_y = flash.YCenter();
      _flash_z = flash.ZCenter();
      _tpc_x   = match.tpc_point.x;
      _tpc_y   = match.tpc_point.y;
      _tpc_z   = match.tpc_point.z;
      _npe     = flash.TotalPE();
      _score   = match.score;
      _flash_time = flash.Time();
      _flash_x = (_flash_time/0.5)*_t2cm; // 0.5 is sampling time
      _flash_y = flash.YCenter();
      
      // if the z-distance is not too large -> assume good match
      // now try and find the michel
      auto matched_flash = _FindOpMichel.FindMichelMatch(*ev_flash,flash);
      std::cout << "matched flash is : " << matched_flash << std::endl;
      if (matched_flash >= 0){
	auto michel_flash = (*ev_flash)[matched_flash];
	_muon_t    = flash.Time();
	_muon_pe   = flash.TotalPE();
	_michel_t  = michel_flash.Time();
	_michel_pe = michel_flash.TotalPE();
	_dz        = fabs(flash.ZCenter() - michel_flash.ZCenter());
	_dy        = fabs(flash.YCenter() - michel_flash.YCenter());
	std::cout << "muon T    = " << _muon_t << std::endl
		  << "michel T  = " << _michel_t << std::endl
		  << "delta T   = " << fabs(_muon_t-_michel_t) << std::endl
		  << "delta Z   = " << _dz << std::endl
		  << "delta Y   = " << _dy << std::endl
		  << "muon PE   = " << _muon_pe << std::endl
		  << "michel PE = " << _michel_pe << std::endl << std::endl;
	_match_tree->Fill();
      }

      _mc_time = _mc_x = _mc_y = _mc_z = -1;
      if(_use_mc) {
	auto const& mct = (*ev_mctrack)[match.tpc_id];
	_mc_time = mct[0].T() * 1.e-3;
	double min_dist = 1e12;
	pt[0] = _tpc_x;
	pt[1] = _tpc_y;
	pt[2] = _tpc_z;
	double min_x=1e9;
	double max_x=0;
	for(size_t i=0; i<mct.size()-1; ++i) {
	  auto const& step1 = mct[i];
	  auto const& step2 = mct[i+1];
	  line.Start(step1.X(),step1.Y(),step1.Z());
	  line.End(step2.X(),step2.Y(),step2.Z());
	  auto const closest_pt = geoalg.ClosestPt(line,pt);
	  double dist = closest_pt.SqDist(pt);
	  if(dist < min_dist) {
	    min_dist = dist;
	    _mc_x = closest_pt[0];
	    _mc_y = closest_pt[1];
	    _mc_z = closest_pt[2];
	  }

	  if(step1.X() < min_x) min_x = step1.X();
	  if(step1.X() > max_x) max_x = step1.X();
	}
	_mc_dx = max_x - min_x;
      }// if we are using mc info
      _tree->Fill();
    }

    return true;
  }

  bool MuonClusterTagger::finalize() {

    if (_fout){
      _fout->cd();
      if (_tree)
	_tree->Write();
      if (_match_tree)
	_match_tree->Write();
    }

    return true;
  }

}
#endif
