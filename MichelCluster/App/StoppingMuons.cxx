#ifndef LARLITE_STOPPINGMUONS_CXX
#define LARLITE_STOPPINGMUONS_CXX

#include "StoppingMuons.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

//stolen from kas

namespace larlite {

  bool StoppingMuons::initialize() {
    
    michel_filter_tree = new TTree("michel_filter_tree","michel_filter_tree");
    michel_filter_tree->Branch("_michel_x",&_michel_x,"michel_x/D");
    michel_filter_tree->Branch("_michel_energy",&_michel_energy,"michel_energy/D");
    michel_filter_tree->Branch("_michel_charge",&_michel_charge,"_michel_charge/D");
    michel_filter_tree->Branch("_michel_det",&_michel_det,"_michel_det/D");
    michel_filter_tree->Branch("_lifetime_correction",&_lifetime_correction,"_lifetime_correction/D");

    event_filter_tree = new TTree("event_filter_tree","Event Filter Tree");
    event_filter_tree->Branch("_run",&_run,"run/I");    
    event_filter_tree->Branch("_subrun",&_subrun,"subrun/I");
    event_filter_tree->Branch("_event",&_event,"event/I");
    event_filter_tree->Branch("_n_michel",&_n_michel,"n_michel/I");
    
    total_evts = 0;
    kept_evts = 0;

    return true;
  }
  
  bool StoppingMuons::analyze(storage_manager* storage) {

    total_evts++;

    // reset nnumber of michels in the event
    _n_michel = 0;

    //Grab the MCShowers
    auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");    
    if(!ev_mctrack) {
      print(larlite::msg::kERROR,__FUNCTION__,Form("Did not find specified data product, mctrack!"));
      return false;
    }  

    _run    = storage->run_id();
    _subrun = storage->subrun_id();
    _event  = storage->event_id();
    
    //If no MCTracks in event, no michels for sure
    if(!ev_mctrack->size()){
      if(_flip) { kept_evts++; return true; }
    }

    bool there_is_michel = false;

    //Loop over MCShowers, ask if they came from a muon decay
    for(auto const& mct : *ev_mctrack){

      if (mct.size() < 2)
	continue;

      //    if ( (mct[mct.size()-1].X() > 0) &&  
      /*

      if( ((mct.MotherPdgCode() == 13) or (mcs.MotherPdgCode() == -13) ) &&
	  ( (mct.Process() == "Decay") ) ) { //or (mcs.Process() == "muMinusCaptureAtRest") ) ) {

	_michel_energy = mcs.Start().E();
	_michel_charge = mcs.Charge(2);
	_michel_det    = mcs.DetProfile().E();
	_michel_x      = mcs.DetProfile().X();
	auto t         =  _michel_x/160.;
	_lifetime_correction = exp(t/3.0);

	michel_filter_tree->Fill();
	there_is_michel = true;
	_n_michel += 1;

      }
      */
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

  bool StoppingMuons::finalize() {

    michel_filter_tree->Write();
    event_filter_tree->Write();
    
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("~~~~StoppingMuonsing~~~~"));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",total_evts,kept_evts));
    
    return true;
  }

  
}
#endif
