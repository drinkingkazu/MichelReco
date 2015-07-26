#ifndef LARLITE_MICHELRECODRIVER_CXX
#define LARLITE_MICHELRECODRIVER_CXX

#include "MichelRecoDriver.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
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

    // Get data products
    auto ev_cluster = storage->get_data<event_cluster>(_producer);

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
			       h.WireID().Wire * 0.3,
			       (h.PeakTime() - 3200) * 0.0802814,
			       hit_index );
    }
    
    // Reaching this point means we have something to process. Prepare.
    _mgr.EventReset();
    
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
				     hit_index );
      }

      // Append a hit-list (cluster) to a manager if not empty
      if(michel_cluster.size())
	//_mgr.Append(michel_cluster);
	_mgr.Append(std::move(michel_cluster));
    }

    _mgr.RegisterAllHits( std::move(all_hits_v) );
    
    // Now process
    _mgr.Process();

    return true;
  }

  bool MichelRecoDriver::finalize() {
    _mgr.Finalize(_fout);
    return true;
  }

}
#endif
