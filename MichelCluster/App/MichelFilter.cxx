#ifndef LARLITE_MICHELFILTER_CXX
#define LARLITE_MICHELFILTER_CXX

#include "MichelFilter.h"
#include "DataFormat/mcshower.h"


//stolen from kas

namespace larlite {

  bool MichelFilter::initialize() {
    
    //n_events_tree = new TTree("n_events_tree","n_events_tree");

    total_evts = 0;
    kept_evts = 0;

    return true;
  }
  
  bool MichelFilter::analyze(storage_manager* storage) {
    total_evts++;
    
    //Grab the MCShowers
    auto ev_mcshower = storage->get_data<event_mcshower>("mcreco");    
    if(!ev_mcshower) {
      print(larlite::msg::kERROR,__FUNCTION__,Form("Did not find specified data product, mcshower!"));
      return false;
    }  
    
    //If no MCShowers in event, no michels for sure
    if(!ev_mcshower->size()){
      return false;
    }
    //Loop over MCShowers, ask if they came from a muon decay
    for(auto const& mcs : *ev_mcshower){
      if((mcs.MotherPdgCode() == 13                &&
	  mcs.Process() == "muMinusCaptureAtRest") &&
	 (mcs.DetProfile().E()/mcs.Start().E()  > 0.5
	  || mcs.DetProfile().E() >= 15)) {
	kept_evts++;
	return true;
      }
    }
    
    //std::cout<<"michel filter returning false!"<<std::endl;
    return false;
  }

  bool MichelFilter::finalize() {

    print(larlite::msg::kNORMAL,__FUNCTION__,Form("~~~~MichelFiltering~~~~"));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    
    return true;
  }

  
}
#endif
