#ifndef LARLITE_TAGSTOPPINGMUON_CXX
#define LARLITE_TAGSTOPPINGMUON_CXX

#include "TagStoppingMuon.h"
#include "LArUtil/GeometryHelper.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "DataFormat/mctrack.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include "MichelCluster/Fmwk/MichelException.h"

namespace larlite {

  TagStoppingMuon::TagStoppingMuon()
    : _algo(nullptr)
  {
    _fout = 0;
    _name = "TagStoppingMuon";
    _clus_producer = "cccluster";
    _minMuonLength = 30;
  }

  bool TagStoppingMuon::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("_tree","Stopping Tree");
    _tree->Branch("_Estop",&_Estop,"Estop/D");
    _tree->Branch("_dS","std::vector<double>",&_dS);
    _tree->Branch("_dQ","std::vector<double>",&_dQ);
    _tree->Branch("_muonMIPdS","std::vector<double>",&_muonMIPdS);
    _tree->Branch("_muonMIPdQ","std::vector<double>",&_muonMIPdQ);
    _tree->Branch("_MIPdS","std::vector<double>",&_MIPdS);
    _tree->Branch("_MIPdQ","std::vector<double>",&_MIPdQ);
    _tree->Branch("_MIPm",&_MIPm,"MIPm/D");
    _tree->Branch("_MIPs",&_MIPs,"MIPs/D");
    _tree->Branch("_MIPendIdx",&_MIPendIdx,"MIPendIdx/I");
    _tree->Branch("_MIPmedian",&_MIPmedian,"MIPmedian/D");
    _tree->Branch("_MIPrms",&_MIPrms,"MIPrms/D");
    _tree->Branch("_MIPrms_corr",&_MIPrms_corr,"MIPrms_corr/D");
    _tree->Branch("_braggIdx",&_braggIdx,"braggIdx/I");
    _tree->Branch("_braggQ",&_braggQ,"braggQ/D");
    _tree->Branch("_braggAmp",&_braggAmp,"braggAmp/D");
    _tree->Branch("_braggArea",&_braggArea,"braggArea/D");
    _tree->Branch("_braggExpected",&_braggExpected,"braggExpected/D");
    _tree->Branch("_braggDistToEnd",&_braggDistToEnd,"braggDistToEnd/D");
    

    return true;
  }
  
  bool TagStoppingMuon::analyze(storage_manager* storage) {

    double w2cm = larutil::GeometryHelper::GetME()->WireToCm();
    double t2cm = larutil::GeometryHelper::GetME()->TimeToCm();

    _dQ.clear();
    _dS.clear();


    // Start by looking at MC info if available
    //Grab the MCTracks
    auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");    
    
    if(!ev_mctrack) {
      print(larlite::msg::kERROR,__FUNCTION__,Form("Did not find specified data product, mctrack!"));
    }  
    //If no MCTracks in event, no muon for sure
    else if(!ev_mctrack->size()){
      print(larlite::msg::kERROR,__FUNCTION__,Form("No mctracks!"));
    }
    else{
      for(auto const& mct : *ev_mctrack){
	if(abs(mct.PdgCode()) == 13){
	  _Estop = mct[mct.size()-1].E()-105.65;
	  break;
	}
      }// for all mctracks
    }
    
    auto ev_cluster = storage->get_data<event_cluster>(_clus_producer);
    
    if (!ev_cluster){
      std::cout << "No cluster data product in this event...exit" << std::endl;
      return false;
    }
    
    // Get hits & association
    event_hit* ev_hit = nullptr;
    auto const& hit_ass_set = storage->find_one_ass(ev_cluster->id(), ev_hit, ev_cluster->name());
    
    // If ev_hit is null, failed to find data or assocaition (shouldn't happen)
    if(!ev_hit || ev_hit->empty()) return false;
    
    // save all hits in the event
    std::vector< ::michel::HitPt > all_hits_v;
    all_hits_v.reserve(ev_hit->size());
    
    for(size_t hit_index=0; hit_index<ev_hit->size(); ++hit_index) {
      auto const& h = (*ev_hit)[hit_index];
      
      all_hits_v.emplace_back( h.Integral(), h.WireID().Wire*w2cm,
			       (h.PeakTime()-3200)*t2cm, h.WireID().Plane);
      
    }// for all hits

    std::vector<michel::MichelCluster> candidateMuons;

    // Loop over clusters & add them to our algorithm manager
    for(auto const& hit_ass : hit_ass_set) {
      
      // Prepare our hit-list representation    
      std::vector< ::michel::HitPt > clusHits;
      clusHits.reserve(hit_ass.size());
      
      // Loop over hits and fill
      for(auto const& hit_index : hit_ass) {
	
	auto const& h = (*ev_hit)[hit_index];

	if(h.WireID().Plane != 2) continue;
	
	clusHits.emplace_back( h.Integral(), h.WireID().Wire * w2cm,
			       (h.PeakTime() - 3200) * t2cm, hit_index,
			       h.WireID().Plane);
      }// for hits in cluster
      
      // make a Cluster object
      if (clusHits.size() <= 3) continue;

      michel::MichelCluster Clus(clusHits,3,3);

      candidateMuons.push_back(Clus);
      
    }// for all clusters

    // select only longest cluster in the event
    size_t muonIdx = michel::kINVALID_SIZE;
    double lenMax = 0;
    for (size_t i=0; i < candidateMuons.size(); i++){
      auto& clus = candidateMuons[i];
      double len = clus._hits[0].SqDist(clus._hits[clus._hits.size()-1]);
      if (len > lenMax){
	lenMax = len;
	muonIdx = i;
      }// if we found the longest muon
    }//for all muons

    if (muonIdx == michel::kINVALID_SIZE) return true;
    
    auto& muon = candidateMuons[muonIdx];

    _algo->ProcessCluster(muon,all_hits_v);

    // get truncated Q vector
    _dQ = muon._t_mean_v;
    // get dS vector
    _dS = muon._s_v;

    //std::cout << "track length: " << _dS[_dS.size()-1] << std::endl;

    // cut on minimum length:
    if (_dS[_dS.size()-1] < _minMuonLength)
      return true;

    auto braggInfo = GetMaxIndex(_dQ);
    // if we find the bragg peak in the 1st half
    // -> flip the order of the vector so that
    // the muon is moving towards the right
    if ( braggInfo.first < (size_t)(_dQ.size()/2.) ){
      std::reverse(_dS.begin(),_dS.end());
      std::reverse(_dQ.begin(),_dQ.end());
      size_t newIdx = _dQ.size()-braggInfo.first-1;
      braggInfo = std::pair<size_t,double>(newIdx,braggInfo.second);
    }
    _braggIdx = braggInfo.first;
    _braggQ   = braggInfo.second;

    _braggDistToEnd = abs( _dS[_braggIdx] - _dS[_dS.size()-1] );

    // get the point that is 15 cm away from the Bragg Peak
    _MIPendIdx = GetMIPendPos(_dS,braggInfo.first,15);
    if (_MIPendIdx == michel::kINVALID_SIZE)
      return true;

    // make a vector of the MIP region only
    _muonMIPdS = std::vector<double>(_dS.begin(),_dS.begin()+_MIPendIdx);
    _muonMIPdQ = std::vector<double>(_dQ.begin(),_dQ.begin()+_MIPendIdx);

    // get the median and RMS of the MIP ADC values
    // make a copy 'cause we don't want the original sorted
    auto muonMIPdQsorted = _muonMIPdQ;
    _MIPmedian = GetMedian(muonMIPdQsorted);
    _MIPrms    = GetRms(_muonMIPdQ);

    // get the hit indices for MIP hits within 1 RMS of median
    std::vector<size_t> MIPindices = GetMIPindices(_muonMIPdQ,_MIPmedian,_MIPrms,1.);
    // get the vectors for these indices
    _MIPdS = GetSubVector(_muonMIPdS,MIPindices);
    _MIPdQ = GetSubVector(_muonMIPdQ,MIPindices);
    // get the linear fit values for these vectors
    auto fit = GetLinearFit(_MIPdS,_MIPdQ);
    _MIPs = fit.first;
    _MIPm = fit.second;

    // corrected RMS (taking into account only truly MIP hits)
    _MIPrms_corr = GetRms(_MIPdS,_MIPdQ,_MIPs,_MIPm);

    // calcualte bragg amplitude/area
    _braggExpected = _MIPm + _MIPs * _braggIdx;
    _braggAmp = _braggQ - _braggExpected;
    _braggArea = GetBraggArea(_dS,_dQ,_MIPendIdx,_braggIdx,_MIPm,_MIPs);
    

    // fill the tree
    _tree->Fill();
    
    return true;
  }

  bool TagStoppingMuon::finalize() {

    if (_fout && _tree)
      _tree->Write();

    return true;
  }

  
  size_t TagStoppingMuon::GetMIPendPos(const std::vector<double>& v,
				       const size_t& max,
				       const double distAsked)
  {

    double distToPeak = v[max];
    
    for (size_t i=0; i < max; i++){
      
      double dist = distToPeak-v[max-i];
      if (dist > distAsked)
	return max-i;
    }

    std::cout << "could not find a point " << distAsked
	      << " away from the bragg peak..." << std::endl;

    return michel::kINVALID_SIZE;
  }

  std::vector<size_t> TagStoppingMuon::GetMIPindices(const std::vector<double>& dQ,
						     const double& median,
						     const double& rms,
						     const double& alpha)
  {

    std::vector<size_t> indices;

    for (size_t n=0; n < dQ.size(); n++){
      if ( (dQ[n] < median+rms) and (dQ[n] > median-rms) )
	indices.push_back(n);
    }

    return indices;
  }


  std::vector<double> TagStoppingMuon::GetSubVector(const std::vector<double> v,
						    const std::vector<size_t> idx)
  {

    std::vector<double> out;

    for (auto& i : idx){
      if (i >= v.size()){
	std::cerr << "\033[93m[ERROR]\033[00m trying to access element larger than vector size..."
		  << std::endl;
	throw michel::MichelException();
      }
      out.push_back(v[i]);
    }
    
    return out;
  }

  double TagStoppingMuon::GetMedian(std::vector<double>& v)
  {

    std::nth_element(v.begin(),v.begin() + v.size()/2, v.end());
    
    return v[v.size()/2];
  }
  

  double TagStoppingMuon::GetAvg(const std::vector<double>& v)
  {

    double avg = 0;
    for (size_t i=0; i < v.size(); i++)
      avg += v[i];

    if (v.size() != 0)
      return avg/((double)(v.size()));
    
    return 0;
  }


  double TagStoppingMuon::GetRms(const std::vector<double>& v,
				 const double avg)
  {

    double rms = 0;
    for (size_t i=0; i < v.size(); i++)
      rms += (v[i]-avg)*(v[i]-avg);

    if (v.size() != 0)
      return sqrt(rms/((double)v.size()-1));
    
    return 0;
  }


  double TagStoppingMuon::GetRms(const std::vector<double>& v)
  {

    double avg = 0;
    for (size_t i=0; i < v.size(); i++)
      avg += v[i];

    if (v.size() != 0)
      avg /= ((double)(v.size()));
    else
      return 0;

    double rms = 0;
    for (size_t i=0; i < v.size(); i++)
      rms += (v[i]-avg)*(v[i]-avg);

    if (v.size() != 0)
      return sqrt(rms/((double)v.size()-1));
    
    return 0;
  }


  double TagStoppingMuon::GetRms(const std::vector<double>& MIPdS,
				 const std::vector<double>& MIPdQ,
				 const double& slope,
				 const double& intercept)
  {

    double rms = 0.;
    
    for (size_t i=0; i < MIPdS.size(); i++)
      rms += (MIPdQ[i] - ( intercept + slope * MIPdS[i] )) * (MIPdQ[i] - ( intercept + slope * MIPdS[i] ));
    
    if (MIPdS.size() == 0)
      return 0;
  
    rms = sqrt(rms/(double)MIPdS.size());
    
    return rms;
  }


  std::pair<size_t,double> TagStoppingMuon::GetMaxIndex(const std::vector<double>& v)
  {

    double max = 0;
    size_t idx = 0;
    for (size_t i=0; i < v.size(); i++)
      if (v[i] > max) { max = v[i]; idx = i; }

    return std::pair<size_t,double>(idx,max);
  }  





  std::pair<double,double> TagStoppingMuon::GetLinearFit(const std::vector<double>& x, 
							 const std::vector<double>& y)
  {

    // make sure both vectors are the same size!
    if (x.size() != y.size()){
      std::cerr << "\033[93m[ERROR]\033[00m dS and dQ vectors do not have the same length"
		<< std::endl;
      throw michel::MichelException();
    }

    // mean of x*y
    double m_xy = 0;
    // mean of x
    double m_x  = 0;
    // mean of y
    double m_y  = 0;
    // mean of x^2
    double m_xx = 0;
    // mean of y^2
    double m_yy = 0;
    
    for (size_t n=0; n < x.size(); n++){

      double xi = x[n];
      double yi = y[n];

      m_xy += xi*yi;
      m_x  += xi;
      m_y  += yi;
      m_xx += xi*xi;
      m_yy += yi*yi;

    }

    double entries = (double)x.size();

    m_xy /= entries;
    m_x  /= entries;
    m_y  /= entries;
    m_xx /= entries;
    m_yy /= entries;

    double slope = (m_xy-m_x*m_y)/(m_xx-m_x*m_x);
    double intercept = m_y - slope * m_x;

    return std::pair<double,double>(slope,intercept);
  }


  double TagStoppingMuon::GetBraggArea(const std::vector<double>& dS,
				       const std::vector<double>& dQ,
				       const size_t& MIPendIdx,
				       const size_t& braggIdx,
				       const double& MIPm,
				       const double& MIPs)
  {

    if (MIPendIdx >= braggIdx){
      std::cerr << "\033[93m[ERROR]\033[00m MIP end index comes after Bragg Peak index...something is wrong..."
		<< std::endl;
      throw michel::MichelException();
    }

    if (braggIdx > dS.size()){
      std::cerr << "\033[93m[ERROR]\033[00m Bragg peak index is larger than vector size...something is wrong..."
		<< std::endl;
      throw michel::MichelException();
    }

    if (dQ.size() != dS.size()){
      std::cerr << "\033[93m[ERROR]\033[00m dS and dQ vectors are not the same length...something is wrong..."
		<< std::endl;
      throw michel::MichelException();
    }

    double area = 0;
    for (size_t i =  MIPendIdx; i < braggIdx; i++)
      area += dQ[i] - ( MIPm + MIPs * dS[i] );

    return area;
  }

}
#endif
