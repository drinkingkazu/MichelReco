#ifndef LARLITE_MICHELMCHITSTUDY_CXX
#define LARLITE_MICHELMCHITSTUDY_CXX

#include "MichelMCHitStudy.h"

#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/simch.h"

namespace larlite {

  MichelMCHitStudy::MichelMCHitStudy()
    :  _mc_tree(nullptr)
    , _mc_hit_tree(nullptr)
  {
    _name="MichelMCHitStudy";
    _fout=0;
    _debug_mcq = false;
    _fill_hit_tree = false;
  }

  bool MichelMCHitStudy::initialize() {

    _mc_hit_tree = new TTree("_mc_hit_tree","MC Hit Information");
    _mc_hit_tree->Branch("_hit_integral",&_hit_integral,"hit_integral/D");
    _mc_hit_tree->Branch("_hit_mc_q",&_hit_mc_q,"hit_mc_q/D");
    _mc_hit_tree->Branch("_event" , &_event  , "event/I");

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
    // Reconstructed energy & Simch stuff
    _mc_tree->Branch("_mc_energy"      , &_mc_energy           , "mc_energy/D" );
    _mc_tree->Branch("_mc_detprof"     , &_mc_detprof          , "mc_detprof/D");
    _mc_tree->Branch("_reco_energy"    , &_reco_energy         , "reco_energy/D");
    // hit by hit simchannel and ADC information [all michel hits]
    _mc_tree->Branch("_michel_hit_fracReco", "std::vector<double>" , &_michel_hit_fracReco);
    _mc_tree->Branch("_michel_hit_QtotReco", "std::vector<double>" , &_michel_hit_QtotReco);
    _mc_tree->Branch("_michel_hit_ADCtotReco", "std::vector<double>" , &_michel_hit_ADCtotReco);
    _mc_tree->Branch("_michel_hit_fracMC", "std::vector<double>" , &_michel_hit_fracMC);
    _mc_tree->Branch("_michel_hit_QtotMC", "std::vector<double>" , &_michel_hit_QtotMC);
    // hit by hit simchannel and ADC information [electron-only hits]
    _mc_tree->Branch("_electron_hit_fracReco", "std::vector<double>" , &_electron_hit_fracReco);
    _mc_tree->Branch("_electron_hit_QtotReco", "std::vector<double>" , &_electron_hit_QtotReco);
    _mc_tree->Branch("_electron_hit_ADCtotReco", "std::vector<double>" , &_electron_hit_ADCtotReco);
    _mc_tree->Branch("_electron_hit_fracMC", "std::vector<double>" , &_electron_hit_fracMC);
    _mc_tree->Branch("_electron_hit_QtotMC", "std::vector<double>" , &_electron_hit_QtotMC);

    _mc_tree->Branch("_QMichelMC", &_QMichelMC, "QMichelMC/D");
    _mc_tree->Branch("_QMichelRecoSimch_all", &_QMichelRecoSimch_all, "QMichelRecoSimch_all/D");
    _mc_tree->Branch("_QMichelRecoSimch_shr", &_QMichelRecoSimch_shr, "QMichelRecoSimch_shr/D");
    _mc_tree->Branch("_QElectronRecoSimch_all", &_QElectronRecoSimch_all, "QElectronRecoSimch_all/D");
    _mc_tree->Branch("_QElectronRecoSimch_shr", &_QElectronRecoSimch_shr, "QElectronRecoSimch_shr/D");
    _mc_tree->Branch("_QMichelShowerMCSimch_all", &_QMichelShowerMCSimch_all, "QMichelShowerMCSimch_all/D");
    _mc_tree->Branch("_QMichelShowerMCSimch_shr", &_QMichelShowerMCSimch_shr, "QMichelShowerMCSimch_shr/D");
    _mc_tree->Branch("_QMichelPartMCSimch_all", &_QMichelPartMCSimch_all, "QMichelPartMCSimch_all/D");
    _mc_tree->Branch("_QMichelPartMCSimch_shr", &_QMichelPartMCSimch_shr, "QMichelPartMCSimch_shr/D");
    _mc_tree->Branch("_f_RecoHitsQ_fromMichelSimch",&_f_RecoHitsQ_fromMichelSimch,"f_RecoHitsQ_fromMichelSimch/D");

    return true;
  }
  
  bool MichelMCHitStudy::analyze(storage_manager* storage) {

    auto ev_mcshower = storage->get_data<event_mcshower>( "mcreco2"  );  
    auto ev_mctrack  = storage->get_data<event_mctrack> ( "mcreco"   );
    auto ev_simch    = storage->get_data<event_simch>   ( "largeant" );
    auto ev_hit      = storage->get_data<event_hit>     ( "gaushit"  );
    auto ev_clus     = storage->get_data<event_cluster> ( "pandoraCosmic" );
    larlite::event_hit *ev_cluster_hit;
    auto const& ass_clus_hit_v = storage->find_one_ass( ev_clus->id(), ev_cluster_hit, ev_clus->name());

    // vector of hit indices to use
    std::vector<unsigned int> hit_idx_v;
    hit_idx_v.push_back(0);

    _event  = storage->event_id();
    _subrun = storage->subrun_id();
    _run    = storage->run_id();

    // are there saved MC showers? if no quit
    if (!ev_mcshower or (ev_mcshower->size() == 0) ){
      std::cout << "No MCShower for module " << _name << " -> quit " << std::endl;
      return false;
    }

    // are there saved MC tracks? if no quit
    if (!ev_mctrack or (ev_mctrack->size() == 0) ){
      std::cout << "No MCTrack for module " << _name << " -> quit " << std::endl;
      return false;
    }

    if (!ev_hit or (ev_hit->size() == 0) ){
      std::cout << "No Hit for module " << _name << " -> quit " << std::endl;
      return false;
    }

    for (size_t i=0; i < ev_mcshower->size(); i++){
      
      auto const& mcs = ev_mcshower->at(i);

      // only keep electron/positron
      if (mcs.PdgCode() == 22)
	continue;

      // only keep > 5 MeV
      if (mcs.Start().E() < 5)
	continue;
      
      // only keep decays
      if ( (mcs.Process() != "Decay") and (mcs.Process() != "muMinusCaptureAtRest") )
	continue;

      Reset();

      // gather up track IDs for all particles in the Michel EM Shower
      std::vector<std::vector<unsigned int> > g4_trackid_v_shr;
      // use only the Michel particle, ignoring radiated photons
      std::vector<std::vector<unsigned int> > g4_trackid_v_part;
      
      _mc_energy  = mcs.Start().E();
      _mc_detprof = mcs.DetProfile().E();
      _QMichelMC = mcs.Charge(2);
      // calculate lifetime correction
      _mc_x  = mcs.Start().X();
      _mc_y  = mcs.Start().Y();
      _mc_z  = mcs.Start().Z();
      _mc_px = mcs.Start().Px();
      _mc_py = mcs.Start().Py();
      _mc_pz = mcs.Start().Pz();
      
      std::vector<unsigned int> id_v;
      id_v.reserve(mcs.DaughterTrackID().size());
      
      for(auto const& id : mcs.DaughterTrackID()) {
	if(id == mcs.TrackID()) continue;
	id_v.push_back(id);
      }
      id_v.push_back(mcs.TrackID());
      g4_trackid_v_shr.push_back(id_v);
      // only for michel particle and no photons
      std::vector<unsigned int> id_v_part = {mcs.TrackID()};
      g4_trackid_v_part.push_back(id_v_part);
      
      if(g4_trackid_v_shr.size() == 0)
	std::cout << "No track IDs found for backtracker !!!" << std::endl;
      
      else{

	try { _BTAlgShower.BuildMap(g4_trackid_v_shr, *ev_simch, *ev_hit, ass_clus_hit_v); }
	catch(...) { std::cout << "\n Exception at building BTAlg map...\n"; }
	try { _BTAlgPart.BuildMap(g4_trackid_v_part, *ev_simch, *ev_hit, ass_clus_hit_v); }
	catch(...) { std::cout << "\n Exception at building BTAlg map...\n"; }

	auto btalgoShower = _BTAlgShower.BTAlg();
	auto btalgoPart   = _BTAlgPart.BTAlg();
	
	std::vector<double> hit_frac_michel;
	std::vector<double> hit_Qtot_michel;
	
	for (size_t hit_index = 0; hit_index < ev_hit->size(); hit_index++){
	  
	  const auto& h = ev_hit->at(hit_index);
	  
	  if(h.WireID().Plane != 2)
	    continue;
	  
	  ::btutil::WireRange_t wire_hit(h.Channel(),h.PeakTime()-4*h.RMS()+2398,h.PeakTime()+4*h.RMS()+2398);
	  
	  // check if this wire-range has already been used
	  auto wire_ranges = getUnUsedWireRange(wire_hit);
	  
	  // loop through all output wire-ranges
	  for (auto const& wire_range : wire_ranges){
	    
	    //if the wire-ranges are empty, continue
	    if (wire_range.start == wire_range.end)
	      continue;

	    if (_debug_mcq)
	      std::cout << "hit wire : " << h.Channel() << " w/ range [" << wire_range.start << ", " << wire_range.end << "]" << " with integral : " << h.Integral() << std::endl;
	    
	    auto partsShower = btalgoShower.MCQ(wire_range);
	    
	    // calculate hit MC information
	    _hit_integral = h.Integral();
	    _hit_mc_q     = partsShower.at(0)+partsShower.at(1);
	    if (_fill_hit_tree)
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
	    
	    auto partsPart   = btalgoPart.MCQ(wire_range);

	    // now add charge from michel original electron only
	    double electron_part = partsPart.at(0);
	    other_part  = partsPart.at(1);
	    hit_frac    = electron_part / ( electron_part + other_part );
	    if (hit_frac > 0.1){
	      _QMichelPartMCSimch_shr += electron_part;
	      _QMichelPartMCSimch_all += electron_part + other_part;
	      _electron_hit_QtotMC.push_back(electron_part);
	      _electron_hit_fracMC.push_back(hit_frac);
	    }

	  }// for all wire-ranges
	}// for all hits
	
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
      
      _mc_tree->Fill();
      
    } // for all reco showers
    
    return true;
  }

  bool MichelMCHitStudy::finalize() {

    _mc_tree->Write();
    _mc_hit_tree->Write();

    return true;
  }

  std::vector<btutil::WireRange_t> MichelMCHitStudy::getUnUsedWireRange(const btutil::WireRange_t& wirerange){

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


  void MichelMCHitStudy::Reset()
  {
    _QMichelRecoSimch_all       = 0.;
    _QMichelRecoSimch_shr       = 0.;
    _QElectronRecoSimch_all     = 0.;
    _QElectronRecoSimch_shr     = 0.;
    _QMichelPartMCSimch_all     = 0.;
    _QMichelPartMCSimch_shr     = 0.;
    _QMichelShowerMCSimch_all   = 0.;
    _QMichelShowerMCSimch_shr   = 0.;
    _michel_hit_fracReco.clear();
    _michel_hit_QtotReco.clear();
    _michel_hit_ADCtotReco.clear();
    _michel_hit_idxReco.clear();
    _michel_hit_fracMC.clear();
    _michel_hit_QtotMC.clear();
    _electron_hit_fracReco.clear();
    _electron_hit_QtotReco.clear();
    _electron_hit_ADCtotReco.clear();
    _electron_hit_idxReco.clear();
    _electron_hit_fracMC.clear();
    _electron_hit_QtotMC.clear();
    _wiremap.clear();

    return;
  }


}
#endif
