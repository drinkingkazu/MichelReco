#ifndef LARLITE_REMOVEMICHEL_CXX
#define LARLITE_REMOVEMICHEL_CXX

#include "RemoveMichel.h"
#include "DataFormat/mcshower.h"


//stolen from kas

namespace larlite {

  bool RemoveMichel::initialize() {
    
    total_evts = 0;
    kept_evts = 0;

    return true;
  }
  
  bool RemoveMichel::analyze(storage_manager* storage) {
    total_evts++;
    
    //Grab the MCShowers
    auto ev_mcshower = storage->get_data<event_mcshower>("mcreco");    
    if(!ev_mcshower) {
      print(larlite::msg::kERROR,__FUNCTION__,Form("Did not find specified data product, mcshower!"));
      kept_evts++;
      return true;
    }  
    
    //If no MCShowers in event, no michels for sure
    if(!ev_mcshower->size()){
      kept_evts++;
      return true;
    }
    
    //Loop over MCShowers, ask if they came from a muon decay
    for(auto const& mcs : *ev_mcshower){
      if(mcs.MotherPdgCode() == 13 &&
	 mcs.Process() == "muMinusCaptureAtRest") {
	//mcs.DetProfile().E()/mcs.Start().E() > 0.5) {
	return false;
      }
    }
    
    kept_evts++;
    return true;
  }
  
  bool RemoveMichel::finalize() {
    
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("~~~~RemoveMichel~~~~"));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));    
    return true;
  }

}
#endif
