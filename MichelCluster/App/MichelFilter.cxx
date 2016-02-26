#ifndef LARLITE_MICHELFILTER_CXX
#define LARLITE_MICHELFILTER_CXX

#include "MichelFilter.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

//stolen from kas

namespace larlite {

  bool MichelFilter::initialize() {
    
    michel_filter_tree = new TTree("michel_filter_tree","michel_filter_tree");
    michel_filter_tree->Branch("_michel_x",&_michel_x,"michel_x/D");
    michel_filter_tree->Branch("_michel_y",&_michel_y,"michel_y/D");
    michel_filter_tree->Branch("_michel_z",&_michel_z,"michel_z/D");
    michel_filter_tree->Branch("_michel_energy",&_michel_energy,"michel_energy/D");
    michel_filter_tree->Branch("_michel_charge",&_michel_charge,"michel_charge/D");
    michel_filter_tree->Branch("_michel_det",&_michel_det,"michel_det/D");
    michel_filter_tree->Branch("_michel_process",&_michel_process);
    michel_filter_tree->Branch("_contained",&_contained,"contained/I");

    event_filter_tree = new TTree("event_filter_tree","Event Filter Tree");
    event_filter_tree->Branch("_run",&_run,"run/I");    
    event_filter_tree->Branch("_subrun",&_subrun,"subrun/I");
    event_filter_tree->Branch("_event",&_event,"event/I");
    event_filter_tree->Branch("_n_michel",&_n_michel,"n_michel/I");
    
    total_evts = 0;
    kept_evts = 0;

    return true;
  }
  
  bool MichelFilter::analyze(storage_manager* storage) {

    total_evts++;

    // reset nnumber of michels in the event
    _n_michel = 0;

    //Grab the MCShowers
    auto ev_mcshower = storage->get_data<event_mcshower>("mcreco");    
    if(!ev_mcshower) {
      print(larlite::msg::kERROR,__FUNCTION__,Form("Did not find specified data product, mcshower!"));
      return false;
    }  

    _run    = storage->run_id();
    _subrun = storage->subrun_id();
    _event  = storage->event_id();
    
    //If no MCShowers in event, no michels for sure
    if(!ev_mcshower->size()){
      if(_flip) { kept_evts++; return true; }
    }

    bool there_is_michel = false;

    //Loop over MCShowers, ask if they came from a muon decay
    for(auto const& mcs : *ev_mcshower){

      if( ((mcs.MotherPdgCode() == 13) or (mcs.MotherPdgCode() == -13) ) &&
	  ( mcs.Process() == "Decay" ) ) {
	//( (mcs.Process() == "Decay") ) ) { //or (mcs.Process() == "muMinusCaptureAtRest") ) ) {
	
	_michel_process   = mcs.Process();
	_michel_energy = mcs.Start().E();
	_michel_charge = mcs.Charge(2);
	_michel_det    = mcs.DetProfile().E();
	_michel_x      = mcs.DetProfile().X();
	_michel_y      = mcs.DetProfile().Y();
	_michel_z      = mcs.DetProfile().Z();
	_contained     = 0;

	// make sure the MCShower is in the detector
	if ( (_michel_x > 0)    && (_michel_x < 256)  &&
	     (_michel_x > -116) && (_michel_y < 116)  &&
	     (_michel_z > 0)    && (_michel_z < 1036) ){
	  _contained = 1;
	  
	  michel_filter_tree->Fill();
	  there_is_michel = true;
	  _n_michel += 1;
	}// if in FV
	  
      }
    } // for all mc showers

    event_filter_tree->Fill();
    
    if (_flip){
      if(_n_michel == 0){
	kept_evts++;
	return true;
      }
    }

    if(_n_michel > 0){
      kept_evts++;
      return true;
    }
    
    return false;
  }

  bool MichelFilter::finalize() {

    michel_filter_tree->Write();
    event_filter_tree->Write();
    
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("~~~~MichelFiltering~~~~"));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    
    return true;
  }

  
}
#endif
