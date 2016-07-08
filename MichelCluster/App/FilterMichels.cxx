#ifndef LARLITE_FILTERMICHELS_CXX
#define LARLITE_FILTERMICHELS_CXX

#include "FilterMichels.h"
#include "DataFormat/cluster.h"

namespace larlite {

  bool FilterMichels::initialize() {

    _totl_events = _pass_events = 0;

    return true;
  }
  
  bool FilterMichels::analyze(storage_manager* storage) {

    _totl_events += 1;

    auto ev_cluster = storage->get_data<event_cluster>(_producer);

    if (!ev_cluster or (ev_cluster->size() == 0) )
      return false;
    
    _pass_events += 1;
    return true;
  }

  bool FilterMichels::finalize() {

    std::cout << std::endl << std::endl;
    std::cout << "***************************************************" << std::endl;
    std::cout << "Searching for clusters w/ producer name : " << _producer << std::endl; 
    std::cout << "Total events scanned : " << _totl_events << std::endl;
    std::cout << "Michel events passed : " << _pass_events << std::endl;
    std::cout << "***************************************************" << std::endl;
    std::cout << std::endl << std::endl;

    return true;
  }

}
#endif
