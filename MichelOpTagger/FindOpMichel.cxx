#ifndef LARLITE_FINDOPMICHEL_CXX
#define LARLITE_FINDOPMICHEL_CXX

#include "FindOpMichel.h"
#include <cmath>

namespace larlite {

  FindOpMichel::FindOpMichel()
    : _tree(nullptr)
  {
    _name = "FindOpMichel";
    _time_window = 5; // usec
    _verbose = false;
  }

  void FindOpMichel::initialize()
  {
    if (_tree) delete _tree;
    _tree = new TTree("_tree","op match tree");
    _tree->Branch("_n_compat",&_n_compat,"n_compat/I");
    _tree->Branch("_max_PE",&_max_PE,"max_PE/D");
    _tree->Branch("_min_d",&_min_d,"min_d/D");
    _tree->Branch("_match_PE",&_match_PE,"match_PE/D");
    _tree->Branch("_match_d",&_match_d,"match_d/D");

    return;
  }

  int FindOpMichel::FindMichelMatch(const std::vector<larlite::opflash>& flashes,
				    const opflash& muon)
  {

    // keep track of the maximum score found
    double max_score = -1;
    int best_match = -1;
    // keep track of the flash with the max num. of PE
    _max_PE = 0;
    // keep track of the number of candidates
    _n_compat = 0;

    // loop through all flashes that need to be searched
    for (size_t n=0; n < flashes.size(); n++){

      // the flash
      auto const& flash = flashes.at(n);

      // skip if the two flashes are identical
      if ( (muon.Time() == flash.Time()) and (muon.TotalPE() == flash.TotalPE()) )
	continue;

      // first check that the two flashes are within the allowed
      // time-window
      if ( (flash.Time() - muon.Time() > _time_window) or (flash.Time() < muon.Time()) )
	continue;

      _n_compat += 1;

      if (flash.TotalPE() > _max_PE)
	_max_PE = flash.TotalPE();

      auto score = MatchScore(muon,flash);
      if (score > max_score) { 
	_min_d = 1./score;
	_match_d = _min_d;
	_match_PE = flash.TotalPE();
	max_score  = score;
	best_match = n;
      }

    }// for all flashes

    _tree->Fill();

    if (_verbose)
      std::cout << "num candidates : " << _n_compat << std::endl;
    
    // return the best match (defaults to -1 if none found)
    return best_match;
  }

  double FindOpMichel::MatchScore(const opflash& muon, const opflash& michel) const
  {

    // return distance in y-z space between flashes
    double dist = sqrt ( ( (muon.YCenter() - michel.YCenter()) * (muon.YCenter() - michel.YCenter()) ) +
			 ( (muon.ZCenter() - michel.ZCenter()) * (muon.ZCenter() - michel.ZCenter()) ) );

    return 1./dist;
  }

}
#endif
