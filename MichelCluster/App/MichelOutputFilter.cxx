#ifndef LARLITE_MICHELOUTPUTFILTER_CXX
#define LARLITE_MICHELOUTPUTFILTER_CXX

#include "MichelOutputFilter.h"
#include "DataFormat/cluster.h"

namespace larlite {

  MichelOutputFilter::MichelOutputFilter()
  {
    _name = "MichelOutputFilter";
    _fout = 0;
    _michel_producer = "michel";
  }

  bool MichelOutputFilter::initialize() {

    return true;
  }
  
  bool MichelOutputFilter::analyze(storage_manager* storage) {

    auto ev_clus = storage->get_data<event_cluster>(_michel_producer);

    if ( (!ev_clus) or (ev_clus->size() == 0) )
      return false;
  
    return true;
  }

  bool MichelOutputFilter::finalize() {

   return true;
  }

}
#endif
