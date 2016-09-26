#ifndef LARLITE_MICHELTRUTHSTUDIES_CXX
#define LARLITE_MICHELTRUTHSTUDIES_CXX

#include "MichelTruthStudies.h"

namespace larlite {

  MichelTruthStudies::MichelTruthStudies()
    : _tree_mc(nullptr)
  {
    _name = "MichelTruthStudies";
    _fout = 0;
  }

  bool MichelTruthStudies::initialize() {

    if (_tree_mc) delete _tree_mc;

    _tree_mc = new TTree("_tree_mc","tree");

    _tree_mc->Branch("_event",&_event,"event/I");
    _tree_mc->Branch("_run",&_run,"run/I");
    _tree_mc->Branch("_subrun",&_subrun,"subrun/I");

    // physics process / PDG / etc information
    _tree_mc->Branch("_pdg",&_pdg,"pdg/I");
    _tree_mc->Branch("_parent_pdg",&_parent_pdg,"parent_pdg/I");
    _tree_mc->Branch("process",&_process);
    _tree_mc->Branch("_parent_end_E",&_parent_end_E,"parent_end_E/D");
    
    _tree_mc->Branch("_mc_X",&_mc_X,"mc_X/D");
    _tree_mc->Branch("_mc_Y",&_mc_Y,"mc_Y/D");
    _tree_mc->Branch("_mc_Z",&_mc_Z,"mc_Z/D");
    _tree_mc->Branch("_mc_T",&_mc_T,"mc_T/D");

    _tree_mc->Branch("_mu_trkID",&_mu_trkID,"mu_trkID/I");

    _tree_mc->Branch("_mc_muon_E",&_mc_muon_E,"mc_muon_E/D");
    _tree_mc->Branch("_mc_muon_px",&_mc_muon_px,"mc_muon_px/D");
    _tree_mc->Branch("_mc_muon_py",&_mc_muon_py,"mc_muon_py/D");
    _tree_mc->Branch("_mc_muon_pz",&_mc_muon_pz,"mc_muon_pz/D");

    _tree_mc->Branch("_mc_michel_E_true",&_mc_michel_E_true,"mc_michel_E_true/D");
    _tree_mc->Branch("_mc_michel_E_edep",&_mc_michel_E_edep,"mc_michel_E_edep/D");
    _tree_mc->Branch("_mc_michel_E_elec",&_mc_michel_E_elec,"mc_michel_E_elec/D");
    _tree_mc->Branch("_mc_michel_px",&_mc_michel_px,"mc_michel_px/D");
    _tree_mc->Branch("_mc_michel_py",&_mc_michel_py,"mc_michel_py/D");
    _tree_mc->Branch("_mc_michel_pz",&_mc_michel_pz,"mc_michel_pz/D");

    return true;
  }
  
  bool MichelTruthStudies::analyze(storage_manager* storage) {

    // load MCShower / MCTracks
    auto ev_mcshower = storage->get_data<event_mcshower>("mcreco2");
    auto ev_mctrack  = storage->get_data<event_mctrack>("mcreco");
    
    _trackIDMap.clear();

    _event  = storage->event_id();
    _subrun = storage->subrun_id();
    _run    = storage->run_id();

    // build ID -> position map for MCTracks
    for (size_t i=0; i < ev_mctrack->size(); i++)
      _trackIDMap[ ev_mctrack->at(i).TrackID() ] = i;
    
    for (size_t i=0; i < ev_mcshower->size(); i++){
      
      auto const& mcsh = ev_mcshower->at(i);
      
      // only keep decays
      if ( (mcsh.Process() != "Decay") and (mcsh.Process() != "muMinusCaptureAtRest") )
	continue;
      
      auto const& e_strt = mcsh.Start();
      
      auto const& muon = ev_mctrack->at( _trackIDMap[ mcsh.MotherTrackID() ] );
      
      _mu_trkID = mcsh.MotherTrackID();
      
      _pdg          = mcsh.PdgCode();
      _parent_pdg   = muon.PdgCode();
      if (muon.size() < 2)
	_parent_end_E = -1;
      else
	_parent_end_E = muon.at( muon.size() - 2).E();
      _process      = mcsh.Process();

      _mc_X = e_strt.X();
      _mc_Y = e_strt.Y();
      _mc_Z = e_strt.Z();
      _mc_T = e_strt.T();
      
      _mc_michel_E_true  = mcsh.Start().E();
      _mc_michel_E_edep  = mcsh.DetProfile().E();
      _mc_michel_E_elec  = mcsh.End().E();
      _mc_michel_px = mcsh.Start().Px();
      _mc_michel_py = mcsh.Start().Py();
      _mc_michel_pz = mcsh.Start().Pz();
      
      double mag = sqrt( (_mc_michel_px * _mc_michel_px) +
			 (_mc_michel_py * _mc_michel_py) +
			 (_mc_michel_pz * _mc_michel_pz) );
      
      _mc_michel_px /= mag;
      _mc_michel_py /= mag;
      _mc_michel_pz /= mag;
      
      _mc_michel_creation_T = mcsh.DetProfile().T();

      _tree_mc->Fill();

    }// for all mcshowers
      
  
    return true;
  }

  bool MichelTruthStudies::finalize() {

    _tree_mc->Write();
  
    return true;
  }

}
#endif
