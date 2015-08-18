#ifndef LARLITE_MICHELFILTER_CXX
#define LARLITE_MICHELFILTER_CXX

#include "MichelFilter.h"
#include "DataFormat/mcshower.h"


//stolen from kas

namespace larlite {

  bool MichelFilter::initialize() {
    
    michel_filter_tree = new TTree("michel_filter_tree","michel_filter_tree");

    michel_filter_tree->Branch("_michel_charge",&_michel_charge,"_michel_charge/D");
    michel_filter_tree->Branch("_michel_det",&_michel_det,"_michel_det/D");
    michel_filter_tree->Branch("_lifetime_correction",&_lifetime_correction,"_lifetime_correction/D");
    
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
      if(_flip) { kept_evts++; return true; }
    }

    bool there_is_michel = false;

    //Loop over MCShowers, ask if they came from a muon decay
    for(auto const& mcs : *ev_mcshower){
      if((mcs.MotherPdgCode() == 13                &&
	  mcs.Process() == "muMinusCaptureAtRest") &&
	 (mcs.DetProfile().E()/mcs.Start().E()  > 0.5
	  || mcs.DetProfile().E() >= 15)) {

	_michel_charge = mcs.Charge(2);
	_michel_det    = mcs.DetProfile().E();
	auto t = mcs.DetProfile().X()/160.0;
	_lifetime_correction = exp(t/3.0);

	michel_filter_tree->Fill();
	there_is_michel = true;
      }
    }
    
    if (_flip){
      if(!there_is_michel) kept_evts++;
      return !there_is_michel;
    }

    if(there_is_michel) kept_evts++;
    return there_is_michel;
  }

  bool MichelFilter::finalize() {

    michel_filter_tree->Write();
    
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("~~~~MichelFiltering~~~~"));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    
    return true;
  }

  
}
#endif
