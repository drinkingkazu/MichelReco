#ifndef LARLITE_STOPPINGFILTER_CXX
#define LARLITE_STOPPINGFILTER_CXX

#include "StoppingFilter.h"
#include "DataFormat/mctrack.h"
#include <cmath>

namespace larlite {

  StoppingFilter::StoppingFilter()
  {
    _name = "StoppingFilter";
    _total_events = 0;
    _kept_events = 0;
    _filter = true;
    _pass   = true;
    _stopping = true;
    _pure = false;
  }

  bool StoppingFilter::initialize() {
    return true;
  }
  
  bool StoppingFilter::analyze(storage_manager* storage) {
    _total_events++;

    _pass = true;

    //Grab the MCShowers
    auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");    
    if(!ev_mctrack) {
      print(larlite::msg::kERROR,__FUNCTION__,Form("Did not find specified data product, mctrack!"));
      _pass = false;
      return false;
    }  
    
    //If no MCTracks in event, no muon for sure
    if(!ev_mctrack->size()){
      _pass = false;
      return false;
    }

    //Loop over MCTracks, ask for muon
    for(auto const& mct : *ev_mctrack){
      if(abs(mct.PdgCode()) == 13){
	if (mct[0].E() > 200){
	  // where is the end point?
	  auto end = mct[mct.size()-1];

	  /*
	  std::cout << "This muon stops at: ["
		    << end.X() << ", " 
		    << end.Y() << ", "
		    << end.Z() << "] with E = "
		    << end.E() << std::endl;
	  */
	  _E = end.E()-105.65;

	  if (!_filter)
	    return true;

	  // if we want a pure sample ignoring
	  // tracks that leave with low energy
	  if (_pure){
	    if ( (_E > 10) and (_E < 100) ){
	      _pass = false;
	      return false;
	    }
	  }

	  bool in = true;
	  if ( (end.X() < 5) or (end.X() > 251) )
	    in = false;
	  if ( (end.Y() < -111) or (end.Y() > 111) )
	    in = false;
	  if ( (end.Z() < 5) or (end.Z() > 1031) )
	    in = false;

	  if (_stopping){
	    if (in == false) _pass = false;
	    return in;
	  }
	  else{
	    if (!in == false) _pass = false;
	    return !in;
	  }
	}// if has 20 MeV dep
      }// if muon
    }// for all tracks 

    _pass = false;
    return false;
  }

  bool StoppingFilter::finalize() {

    print(larlite::msg::kNORMAL,__FUNCTION__,Form("~~~~StoppingFiltering~~~~"));
    print(larlite::msg::kNORMAL,__FUNCTION__,Form("Total events considered = %zu, kept events = %zu.",_total_events,_kept_events));
    return true;
  }

  
}
#endif
