#ifndef LARLITE_MICHELMCSTUDY_CXX
#define LARLITE_MICHELMCSTUDY_CXX

#include "MichelMCStudy.h"

#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/simch.h"

namespace larlite {

  MichelMCStudy::MichelMCStudy()
    :  _mc_tree(nullptr)
    , _mc_hit_tree(nullptr)
  {
    _name="MichelMCStudy";
    _fout=0;
  }

  bool MichelMCStudy::initialize() {

    _mc_hit_tree = new TTree("_mc_hit_tree","MC Hit Information");
    _mc_hit_tree->Branch("_hit_integral",&_hit_integral,"hit_integral/D");
    _mc_hit_tree->Branch("_hit_mc_q",&_hit_mc_q,"hit_mc_q/D");

    // MC Michel & Muon Information (& Simchannel stuff)
    _mc_tree = new TTree("_mc_tree","MC comparison TTree");
    // event information
    _mc_tree->Branch("_run"   , &_run    , "run/I");
    _mc_tree->Branch("_subrun", &_subrun , "subrun/I");
    _mc_tree->Branch("_event" , &_event  , "event/I");
    // michel MC information
    _mc_tree->Branch("_mc_x",  &_mc_x,  "mc_x/D" );
    _mc_tree->Branch("_mc_y",  &_mc_y,  "mc_y/D" );
    _mc_tree->Branch("_mc_z",  &_mc_z,  "mc_z/D" );
    _mc_tree->Branch("_mc_px", &_mc_px, "mc_px/D");
    _mc_tree->Branch("_mc_py", &_mc_py, "mc_py/D");
    _mc_tree->Branch("_mc_pz", &_mc_pz, "mc_pz/D");
    _mc_tree->Branch("_mc_yplane_angle",&_mc_yplane_angle,"mc_yplane_angle/D");
    // muon MC information
    _mc_tree->Branch("_mu_x",  &_mu_x,  "mu_x/D" );
    _mc_tree->Branch("_mu_y",  &_mu_y,  "mu_y/D" );
    _mc_tree->Branch("_mu_z",  &_mu_z,  "mu_z/D" );
    _mc_tree->Branch("_mu_px", &_mu_px, "mu_px/D");
    _mc_tree->Branch("_mu_py", &_mu_py, "mu_py/D");
    _mc_tree->Branch("_mu_pz", &_mu_pz, "mu_pz/D");
    _mc_tree->Branch("_mu_yplane_angle",&_mu_yplane_angle,"mu_yplane_angle/D");
    // 3D dot product between michel and muon direction
    _mc_tree->Branch("_mu_michel_3dangle",&_mu_michel_3dangle,"mu_michel_3dangle/D");
    // Reconstructed energy & Simch stuff
    _mc_tree->Branch("_mc_energy"      , &_mc_energy           , "mc_energy/D");
    _mc_tree->Branch("_reco_energy"    , &_reco_energy         , "reco_energy/D");
    _mc_tree->Branch("_michel_hit_fracReco", "std::vector<double>" , &_michel_hit_fracReco);
    _mc_tree->Branch("_michel_hit_QtotReco", "std::vector<double>" , &_michel_hit_QtotReco);
    _mc_tree->Branch("_michel_hit_idxReco",  "std::vector<double>" , &_michel_hit_idxReco );
    _mc_tree->Branch("_michel_hit_fracMC", "std::vector<double>" , &_michel_hit_fracMC);
    _mc_tree->Branch("_michel_hit_QtotMC", "std::vector<double>" , &_michel_hit_QtotMC);
    _mc_tree->Branch("_QMichelMC", &_QMichelMC, "QMichelMC/D");
    _mc_tree->Branch("_QMichelRecoSimch_all", &_QMichelRecoSimch_all, "QMichelRecoSimch_all/D");
    _mc_tree->Branch("_QMichelRecoSimch_shr", &_QMichelRecoSimch_shr, "QMichelRecoSimch_shr/D");
    _mc_tree->Branch("_QMichelShowerMCSimch_all", &_QMichelShowerMCSimch_all, "QMichelShowerMCSimch_all/D");
    _mc_tree->Branch("_QMichelShowerMCSimch_shr", &_QMichelShowerMCSimch_shr, "QMichelShowerMCSimch_shr/D");
    _mc_tree->Branch("_QMichelPartMCSimch_all", &_QMichelPartMCSimch_all, "QMichelPartMCSimch_all/D");
    _mc_tree->Branch("_QMichelPartMCSimch_shr", &_QMichelPartMCSimch_shr, "QMichelPartMCSimch_shr/D");
    _mc_tree->Branch("_f_RecoHitsQ_fromMichelSimch",&_f_RecoHitsQ_fromMichelSimch,"f_RecoHitsQ_fromMichelSimch/D");

    return true;
  }
  
  bool MichelMCStudy::analyze(storage_manager* storage) {

    auto ev_michel   = storage->get_data<event_cluster> ( "michel"   );  
    auto ev_mcshower = storage->get_data<event_mcshower>( "mcreco"   );  
    auto ev_mctrack  = storage->get_data<event_mctrack> ( "mcreco"   );
    auto ev_simch    = storage->get_data<event_simch>   ( "largeant" );

    // get hits associated with michel clusters
    larlite::event_hit *ev_hit = nullptr;
    auto const& ass_clus_hit_v = storage->find_one_ass(ev_michel->id(), ev_hit, ev_michel->name());

    _QMichelRecoSimch_all     = 0.;
    _QMichelRecoSimch_shr     = 0.;
    _QMichelPartMCSimch_all   = 0.;
    _QMichelPartMCSimch_shr   = 0.;
    _QMichelShowerMCSimch_all = 0.;
    _QMichelShowerMCSimch_shr = 0.;
    _michel_hit_fracReco.clear();
    _michel_hit_QtotReco.clear();
    _michel_hit_idxReco.clear();
    _michel_hit_fracMC.clear();
    _michel_hit_QtotMC.clear();

    _wiremap.clear();

    /*

    // get the 
    
    bool shower_exists = false;
    
    _reco_energy = 0;
    _mc_energy = -1;
    
    //
    // Get the MeV scale energy for this shower
    
    for(auto const& mcs : *ev_mcshower){
      
      if ( ( ( mcs.MotherPdgCode() == 13 ) || ( mcs.MotherPdgCode() == -13 ) ) &&              // parent is a muon
	   ( ( mcs.Process() == "muMinusCaptureAtRest" ) || ( mcs.Process() == "Decay" ) ) &&  // process is consistent with michel decay
	   ( mcs.DetProfile().E()/mcs.Start().E()  > 0.5 ) ) {                                 // at least 50% of energy is contained
	
	shower_exists = true;
      }
    }
    
    //
    // If there is an MC shower in this event, backtrack reconstructed
    // hits and find charge deposition
    
    bool reco_michel_exists = true;
    if (!ev_michel or (ev_michel->size() == 0) )
      reco_michel_exists = false;
    
    if(shower_exists){// && reco_michel_exists ) { 
      
      // You will probably complain about this but I need to use backtracker
      // how well we can self contain the class in current repo is another story
      
      auto mc_energy_min = 0;  // MeV 
      auto mc_energy_max = 65; // MeV 
      
      // gather up track IDs for all particles in the Michel EM Shower
      std::vector<std::vector<unsigned int> > g4_trackid_v_shr;
      // use only the Michel particle, ignoring radiated photons
      std::vector<std::vector<unsigned int> > g4_trackid_v_part;
      
      for(size_t mc_index=0; mc_index < ev_mcshower->size(); ++mc_index) {
	auto const& mcs = (*ev_mcshower)[mc_index];
	
	
	if ( ( ( mcs.MotherPdgCode() == 13 ) || ( mcs.MotherPdgCode() == -13 ) ) &&              // parent is a muon
	     ( ( mcs.Process() == "muMinusCaptureAtRest" ) || ( mcs.Process() == "Decay" ) ) &&  // process is consistent with michel decay
	     ( mcs.DetProfile().E()/mcs.Start().E()  > 0.5 ) ) {                                 // at least 50% of energy is contained
	  
	  _mc_energy = mcs.DetProfile().E();
	  _QMichelMC = mcs.Charge(2);
	  // calculate lifetime correction
	  _mc_x  = mcs.Start().X();
	  _mc_y  = mcs.Start().Y();
	  _mc_z  = mcs.Start().Z();
	  _mc_px = mcs.Start().Px();
	  _mc_py = mcs.Start().Py();
	  _mc_pz = mcs.Start().Pz();
	  double mag = sqrt( (_mc_px*_mc_px) + (_mc_py*_mc_py) + (_mc_pz*_mc_pz) );	    
	  _mc_px /= mag;
	  _mc_py /= mag;
	  _mc_pz /= mag;
	  // 2D angle on Y plane view (in x-z space)
	  _mc_yplane_angle = _mc_pz / sqrt ( _mc_px*_mc_px + _mc_pz*_mc_pz );
	  
	  // find the MCTrack that produced this Michel electron
	  for (auto const& mct : *ev_mctrack){
	    if ( (mcs.MotherTrackID() == mct.TrackID()) and ( abs(mct.PdgCode()) == 13 ) ){
	      // if only 1 step: continue
	      if (mct.size() < 2)
		continue;
	      _mu_energy  = mct[0].E();
	      _mu_px = mct[0].Px();
	      _mu_py = mct[0].Py();
	      _mu_pz = mct[0].Pz();
	      mag = sqrt( (_mu_px*_mu_px) + (_mu_py*_mu_py) + (_mu_pz*_mu_pz) );
	      _mu_px /= mag;
	      _mu_py /= mag;
	      _mu_pz /= mag;
	      // 2D angle on Y plane view (in x-z space)
	      _mu_yplane_angle = _mu_pz / sqrt ( _mu_px*_mu_px + _mu_pz*_mu_pz );
	      // muon - michel 3D angle
	      _mu_michel_3dangle = _mu_px*_mc_px + _mu_py*_mc_py + _mu_pz*_mc_pz;
	    }// if we found the parent muon
	  }// for all mctracks
	  
	  double energy = mcs.DetProfile().E();
	  
	  std::vector<unsigned int> id_v;
	  id_v.reserve(mcs.DaughterTrackID().size());
	  
	  if( mc_energy_min < energy && energy < mc_energy_max ) {
	    for(auto const& id : mcs.DaughterTrackID()) {
	      if(id == mcs.TrackID()) continue;
	      id_v.push_back(id);
	    }
	    id_v.push_back(mcs.TrackID());
	    g4_trackid_v_shr.push_back(id_v);
	    // only for michel particle and no photons
	    std::vector<unsigned int> id_v_part = {mcs.TrackID()};
	    g4_trackid_v_part.push_back(id_v_part);
	  }
	}
      }
      
      if(g4_trackid_v_shr.size() == 0)
	std::cout << "No track IDs found for backtracker !!!" << std::endl;
      
      else{
	
	try { _BTAlgShower.BuildMap(g4_trackid_v_shr, *ev_simch, *ev_hit, hit_ass_set); }
	catch(...) { std::cout << "\n Exception at building BTAlg map...\n"; }
	try { _BTAlgPart.BuildMap(g4_trackid_v_part, *ev_simch, *ev_hit, hit_ass_set); }
	catch(...) { std::cout << "\n Exception at building BTAlg map...\n"; }
	
	auto btalgoShower = _BTAlgShower.BTAlg();
	auto btalgoPart   = _BTAlgPart.BTAlg();
	
	std::vector<double> hit_frac_michel;
	std::vector<double> hit_Qtot_michel;
	
	// hits used for michel
	int michel_hits = 0;
	
	for (size_t hit_index = 0; hit_index < ev_hit->size(); hit_index++){
	  
	  const auto& h = ev_hit->at(hit_index);
	  
	  if(h.WireID().Plane != 2)
	    continue;
	  
	  //::btutil::WireRange_t wire_hit(h.Channel(),h.StartTick()+3050,h.EndTick()+3050);
	  ::btutil::WireRange_t wire_hit(h.Channel(),h.PeakTime()-3*h.RMS()+3050,h.PeakTime()+3*h.RMS()+3050);
	  
	  // check if this wire-range has already been used
	  auto wire_ranges = getUnUsedWireRange(wire_hit);
	  
	  // loop through all output wire-ranges
	  for (auto const& wire_range : wire_ranges){
	    
	    //if the wire-ranges are empty, continue
	    if (wire_range.start == wire_range.end)
	      continue;
	    
	    if (_debug_mcq)
	      std::cout << "hit wire : " << h.Channel() << " w/ range [" << wire_range.start << ", " << wire_range.end << "]" << " with integral : " << h.Integral() << std::endl;
	    
	    // offending malloc if _BTAlg is not class member...
	    auto partsShower = btalgoShower.MCQ(wire_range);
	    auto partsPart   = btalgoPart.MCQ(wire_range);
	    
	    // calculate hit MC information
	    _hit_integral = h.Integral();
	    _hit_mc_q     = partsShower.at(0)+partsShower.at(1);
	    _mc_hit_tree->Fill();
	    
	    
	    double michel_part = partsShower.at(0);
	    double other_part  = partsShower.at(1);
	    double hit_frac    = michel_part / ( michel_part + other_part );
	    if (_debug_mcq)
	      std::cout << "hit " << hit_index << " has " << michel_part+other_part << " and " << hit_frac << " michel contribution..." << std::endl;
	    if (hit_frac > 0.1){
	      _QMichelShowerMCSimch_shr += michel_part;
	      _QMichelShowerMCSimch_all += (michel_part+other_part);
	      _michel_hit_QtotMC.push_back(michel_part);
	      _michel_hit_fracMC.push_back(hit_frac);
	    }
	    
	    // if no michel is reconstructed, skip this part
	    if (reco_michel_exists){
	      // get the reconstructed michel electron
	      if (ev_michel and (ev_michel->size() == 1) ){
		auto michel_cluster = ev_michel->at(0);
		// get list of hits associated to the michel cluster
		auto michel_hit_ass_v = ass_michel_hit_v[0];
		// this list has the indices of the hits associated with the michel
		for (aut const& h_idx : michel_hit_ass_v){
		  if (h_idx == hit_index){
		    if (_debug_mcq)
		      std::cout << "...this hit is michel" << std::endl;
		    michel_hits += 1;
		    auto const& hit = ev_hit->at(h_idx);
		    _reco_charge += hit.Integral();
		    _QMichelRecoSimch_all += (michel_part+other_part);
		    _QMichelRecoSimch_shr += michel_part;
		    _michel_hit_QtotReco.push_back(michel_part);
		    _michel_hit_fracReco.push_back(hit_frac);
		    _michel_hit_idxReco.push_back(h._id);
		  }// if the index maches
	      }// for all michel hits
	      }// if there is a michel cluster reconstructed
	    }// if a reconstructed michel exists
	    
	    // now add charge from michel original electron only
	    michel_part = partsPart.at(0);
	    other_part  = partsPart.at(1);
	    hit_frac    = michel_part / ( michel_part + other_part );
	    if (hit_frac > 0.1){
	      _QMichelPartMCSimch_shr += michel_part;
	      _QMichelPartMCSimch_all += michel_part + other_part;
	    }
	    
	  } // for all wire-ranges
	} // for all hits
	
	_f_RecoHitsQ_fromMichelSimch = _QMichelRecoSimch_shr / _QMichelRecoSimch_all;
	
	if (_debug_mcq){
	  std::cout << std::endl
		    << "QMichel Reconstructed (all) (simch) : " << _QMichelRecoSimch_all << std::endl
		    << "QMichel Reconstructed (shr) (simch) : " << _QMichelRecoSimch_shr << std::endl
		    << "QMichel Shower MC (all) (simch)     : " << _QMichelShowerMCSimch_all << std::endl
		    << "QMichel Shower MC (shr) (simch)     : " << _QMichelShowerMCSimch_shr << std::endl
		    << "QMichel Part MC (all) (simch)       : " << _QMichelPartMCSimch_all << std::endl
		    << "QMichel Part MC (shr) (simch)       : " << _QMichelPartMCSimch_shr << std::endl
		    << "Frac of Q actually from the michel  : " << _f_RecoHitsQ_fromMichelSimch << std::endl
		    << "Frac of  MC Shower reco'd in hits   : " << _QMichelRecoSimch_all/_QMichelShowerMCSimch_shr << std::endl
		    << "Frac of  MC Part   reco'd in hits   : " << _QMichelRecoSimch_all/_QMichelPartMCSimch_shr << std::endl
		    << std::endl << std::endl;
	}
	
      }// if trackIDs are found for the back-tracker
      // How do we really know one of our hits is michel, well it has higher fraction of
      // number of electrons from michel than "not", go ahead and swap for writeout, if this vector
      // is empty in TTree then there was no MC shower
      
      
    } // if we found an mcshower
    
    _mc_tree->Fill();

    */

    return true;
  }

  bool MichelMCStudy::finalize() {

    _mc_tree->Write();
    _mc_hit_tree->Write();

    return true;
  }


  std::vector<btutil::WireRange_t> MichelMCStudy::getUnUsedWireRange(const btutil::WireRange_t& wirerange){

    //std::cout << "searching for already used wire-range @ chan " << wirerange.ch << " in range [" << wirerange.start << ", " << wirerange.end << "]" << std::endl;

    // return a vector of WireRanges
    std::vector<btutil::WireRange_t> out_ranges;

    // is this wire in the map?
    // YES
    if ( _wiremap.find( wirerange.ch ) == _wiremap.end() ){
      // no -> add to the map but don't edit anything
      //std::cout << "\t this channel has not been added..." << std::endl;
      std::pair<double,double> tickrange = std::make_pair(wirerange.start,wirerange.end);
      std::vector< std::pair<double,double> > tickranges = {tickrange};
      _wiremap[ wirerange.ch ] = tickranges;
      // return the same exact wire-range
      //std::cout << "\t adding wire-range @ chan "<< wirerange.ch << " in range [" << wirerange.start << ", " << wirerange.end << "]" << std::endl;
      out_ranges.push_back( wirerange );
    }
    // NO
    else{
      // the channel is there, but is the same time-range?
      auto& wire_ranges = _wiremap[ wirerange.ch ];
      // find if these ticks are included
      // has a wire-range overlap been found?
      bool found = false;
      for (size_t i=0; i < wire_ranges.size(); i++){
	auto wire_range = wire_ranges[i];
	// if wirerange is before 
	if ( wirerange.end < wire_range.first )
	  continue;
	// if wirarange is after
	if ( wirerange.start > wire_range.second )
	  continue;
	// if we made it this far -> there is overlap
	// find the overlap region
	// the two completely overlap, return an empty wire-range
	if ( (wirerange.start >= wire_range.first) and (wirerange.end <= wire_range.second) ){
	  found = true;
	  btutil::WireRange_t newrange( wirerange.ch, 0, 0);
	  //std::cout << "\t adding wire-range @ chan "<< newrange.ch << " in range [" << newrange.start << ", " << newrange.end << "]" << std::endl;
	  out_ranges.push_back ( newrange );
	}
	// if we overflow on the low-end
	if ( wirerange.start < wire_range.first ){
	  found = true;
	  btutil::WireRange_t newrange( wirerange.ch, wirerange.start, wire_range.first);
	  //std::cout << "\t adding wire-range @ chan "<< newrange.ch << " in range [" << newrange.start << ", " << newrange.end << "]" << std::endl;
	  out_ranges.push_back (newrange);
	  wire_ranges[i].first = wirerange.start;
	}
	// if we overlfow on the high-end
	if ( wirerange.end > wire_range.second ){
	  found = true;
	  btutil::WireRange_t newrange( wirerange.ch, wire_range.second, wirerange.end);
	  //std::cout << "\t adding wire-range @ chan "<< newrange.ch << " in range [" << newrange.start << ", " << newrange.end << "]" << std::endl;
	  out_ranges.push_back (newrange);
	  wire_ranges[i].second = wirerange.end;
	}
      }// for all wire-ranges already-found on this channel
      if (!found){
	// return the same exact wire-rang
	//std::cout << "\t this wire-range has not been added..." << std::endl;
	//std::cout << "\t adding wire-range @ chan "<< wirerange.ch << " in range [" << wirerange.start << ", " << wirerange.end << "]" << std::endl;
	out_ranges.push_back( wirerange );
	std::pair<double,double> tickrange = std::make_pair(wirerange.start,wirerange.end);
	std::vector< std::pair<double,double> > tickranges = {tickrange};
	_wiremap[ wirerange.ch ] = tickranges;
      }
    }// if there is some overlap
    
    return out_ranges;
  }


}
#endif
