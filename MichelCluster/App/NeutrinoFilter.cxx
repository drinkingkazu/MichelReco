#ifndef LARLITE_NEUTRINOFILTER_CXX
#define LARLITE_NEUTRINOFILTER_CXX

#include "NeutrinoFilter.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool NeutrinoFilter::initialize() {
    
    _tot = 0;
    _pass = 0;

    return true;
  }
  
  bool NeutrinoFilter::analyze(storage_manager* storage) {

    _tot += 1;

    //Grab the MCTruth
    auto ev_mctruth = storage->get_data<event_mctruth>("generator");
    if (!ev_mctruth) {
      print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
      return false;
    }
  
    // Require exactly one neutrino interaction
    if (ev_mctruth->size() != 1) {
      print(larlite::msg::kINFO, __FUNCTION__, Form("ev_mctruth size is not 1!"));
      return false;
    }

    // get neutrino interaction position:
    auto const& nupos = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position();

    bool contained = ( (nupos[0] > -10)  && (nupos[0] < 266) && 
		       (nupos[1] > -126) && (nupos[1] < 126) &&
		       (nupos[2] > -10)  && (nupos[2] < 1040) );
    
    if (contained)
      return false;

    _pass += 1;
      
    return true;
  }

  bool NeutrinoFilter::finalize() {

    std::cout << "FILTERING" << std::endl << std::endl
	      << "Passed       : " << _pass << std::endl
	      << "PASSING FRAC : " << ((float)_pass) / _tot
	      << std::endl << std::endl;

    return true;
  }

}
#endif
