#ifndef LARLITE_FINDOPMICHEL_CXX
#define LARLITE_FINDOPMICHEL_CXX

#include "FindOpMichel.h"
#include <cmath>

namespace larlite {

  FindOpMichel::FindOpMichel()
  {
    _name = "FindOpMichel";
    _time_window = 5; // usec
  }

  int FindOpMichel::FindMichelMatch(const std::vector<larlite::opflash>& flashes,
				    const opflash& muon) const
  {

    // keep track of the maximum score found
    double max_score = -1;
    int best_match = -1;

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

      auto score = MatchScore(muon,flash);
      if (score > max_score) { 
	max_score = score;
	best_match =n;
      }

    }// for all flashes
    
    // return the best match (defaults to -1 if none found)
    return best_match;
  }

  double FindOpMichel::MatchScore(const opflash& muon, const opflash& michel) const
  {

    // return distance in y-z space between flashes
    double dist = sqrt ( ( (muon.YCenter() - michel.YCenter()) * (muon.YCenter() - michel.YCenter()) ) +
			 ( (muon.ZCenter() - michel.ZCenter()) * (muon.ZCenter() - michel.ZCenter()) ) );

    return dist;
  }

}
#endif
