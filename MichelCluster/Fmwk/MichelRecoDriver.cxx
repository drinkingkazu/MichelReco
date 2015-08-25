#ifndef LARLITE_MICHELRECODRIVER_CXX
#define LARLITE_MICHELRECODRIVER_CXX

#include "MichelRecoDriver.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "LArUtil/GeometryUtilities.h"
#include "LArUtil/LArProperties.h"

namespace larlite {

  bool MichelRecoDriver::initialize() {

    if(_producer.empty()) {
      print(msg::kERROR,__FUNCTION__,"Input data product producer name not specified!");
      throw std::exception();
    }
    
    _tree = new TTree("tree","");

    _mgr.Initialize();

    return true;
  }
  
  bool MichelRecoDriver::analyze(storage_manager* storage) {

    // use instances of LArUtil and GeometryUtilities
    // for (w,t) -> (cm, cm) conversion
    // wire->cm
    double w2cm = larutil::GeometryUtilities::GetME()->WireToCm();
    // time->cm (accounting for different operating voltages)
    double driftVel = larutil::LArProperties::GetME()->DriftVelocity(_Efield,87); // [cm/us]
    // tick width in time
   double tickWidth = 0.5; // [us]
    double t2cm = tickWidth*driftVel;


    // Get data products
    auto ev_cluster = storage->get_data<event_cluster>(_producer);

    // get ID information
    auto run = storage->get_data<event_cluster>(_producer)->run();
    auto subrun = storage->get_data<event_cluster>(_producer)->subrun();
    auto event = storage->get_data<event_cluster>(_producer)->event_id();
    
    michel::EventID id;
    id.run = run;
    id.subrun = subrun;
    id.event = event;

    if(!ev_cluster || ev_cluster->empty()) return false;

    // Get hits & association
    event_hit* ev_hit = nullptr;
    auto const& hit_ass_set = storage->find_one_ass(ev_cluster->id(), ev_hit, ev_cluster->name());

    // If ev_hit is null, failed to find data or assocaition (shouldn't happen)
    if(!ev_hit || ev_hit->empty()) return false;

    // Tracker for used-hit index
    std::vector< ::michel::HitPt > all_hits_v;
    all_hits_v.reserve(ev_hit->size());
    for(size_t hit_index=0; hit_index<ev_hit->size(); ++hit_index) {
      auto const& h = (*ev_hit)[hit_index];
      all_hits_v.emplace_back( h.Integral(),
			       h.WireID().Wire * w2cm,
			       (h.PeakTime() - 3200) * t2cm,
			       hit_index,
			       h.WireID().Plane);
    }
    
    // Reaching this point means we have something to process. Prepare.
    _mgr.EventReset();

    /// set event info
    _mgr.SetEventInfo(id);
        
    // Loop over clusters & add them to our algorithm manager
    for(auto const& hit_ass : hit_ass_set) {

      // Prepare our hit-list representation    
      std::vector< ::michel::HitPt > michel_cluster;
      michel_cluster.reserve(hit_ass.size());

      // Loop over hits and fill
      for(auto const& hit_index : hit_ass) {

	auto const& h = (*ev_hit)[hit_index];
	
	if(h.WireID().Plane != 2) continue;
	
	michel_cluster.emplace_back( h.Integral(),
				     h.WireID().Wire * 0.3,
				     (h.PeakTime() - 3200) * 0.0802814,
				     hit_index ,
				     h.WireID().Plane);
      }
      
      // Append a hit-list (cluster) to a manager if not empty
      if(michel_cluster.size())
	_mgr.Append(std::move(michel_cluster));
    }

    _mgr.RegisterAllHits( std::move(all_hits_v) );
    
    // Now process
    _mgr.Process();

    // We may want to save the michels we found as clusters
    // if this option is selected, take the michel hits
    // and group them in "michel" clusters
    if (_save_clusters){

      auto michel_cluster = storage->get_data<event_cluster>("michel");
      // create association object for clusters
      auto cluster_ass_v = storage->get_data<event_ass>(michel_cluster->name());
      // set event ID through storage manager
      storage->set_id(storage->get_data<event_cluster>(_producer)->run(),
		      storage->get_data<event_cluster>(_producer)->subrun(),
		      storage->get_data<event_cluster>(_producer)->event_id()); 

      // create a vector where to store the cluster -> hit association
      std::vector<std::vector<unsigned int> > clus_hit_ass_v;
      // we need to loop over the michel hits and add them to new clusters

      // get the vector of MichelCluster objects
      auto const& michels = _mgr.GetResult();
      // for each get the hits associated
      for (auto const& michelClus : michels){
	// prepare an empty cluster
	larlite::cluster clus;
	michel_cluster->push_back(clus);
	// get the hits (in "michel" notation) for this cluster
	auto const& michel = michelClus._michel; // this is a vector of HitPt
	// michel is a list of HitPt
	// each one has a HitID_t unique ID which we can use
	// to trace it back to the larlite hit it came from
	// the hit index is indeed the position inside the larlite::hit vector
	// we basically need to create an association between the cluster
	// and all the hits in this michel
	// empty vector where to store hits associated for this specific cluster
	std::vector<unsigned int> clus_hits;
	for (auto const& michel_hit : michel)
	  clus_hits.push_back(michel_hit._id);
	// add this vector to the associations
	if (clus_hits.size() > 3)
	  clus_hit_ass_v.push_back(clus_hits);
	
      }// loop over all found michels
      // now save the association information
      cluster_ass_v->set_association(michel_cluster->id(),product_id(data::kHit,ev_hit->name()),clus_hit_ass_v);
    }// if we should save cluster

    return true;
  }

  bool MichelRecoDriver::finalize() {
    _mgr.Finalize(_fout);
    return true;
  }

}
#endif
