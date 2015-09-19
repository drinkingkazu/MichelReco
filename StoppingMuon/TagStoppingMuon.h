/**
 * \file TagStoppingMuon.h
 *
 * \ingroup StoppingMuon
 * 
 * \brief Class def header for a class TagStoppingMuon
 *
 * @author david caratelli
 */

/** \addtogroup StoppingMuon

    @{*/

#ifndef LARLITE_TAGSTOPPINGMUON_H
#define LARLITE_TAGSTOPPINGMUON_H

#include "Analysis/ana_base.h"
#include "MichelCluster/Fmwk/MichelCluster.h"
#include "MichelCluster/Fmwk/BaseMichelAlgo.h"
#include "MichelCluster/Fmwk/MichelConstants.h"
#include "TTree.h"

namespace larlite {
  /**
     \class TagStoppingMuon
     User custom analysis class made by SHELL_USER_NAME
   */
  class TagStoppingMuon : public ana_base{
  
  public:

    /// Default constructor
    TagStoppingMuon();

    /// Default destructor
    virtual ~TagStoppingMuon(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetClusterProducer(std::string prod) { _clus_producer = prod; }

    void SetAlgo(michel::BaseMichelAlgo* algo) { _algo = algo; }

    // getter for truncated charge vector
    std::vector<double> GetTruncQ() { return _dQ; }
    // getter for dS vector
    std::vector<double> GetdS() { return _dS; }

    // set minimum muon length
    void SetMinMuonLength(double d) { _minMuonLength = d; }

  protected:

    std::string _clus_producer;
    
    michel::BaseMichelAlgo* _algo;



    // settable parameter : minimum muon length required
    double _minMuonLength;

    // get the position marked as the end of the MIP region
    // this is 15 cm before the bragg peak is reached
    size_t GetMIPendPos(const std::vector<double>& v,
			const size_t& max,
			const double distAsked);

    // get indices for hits considered MIP (within some alpha*RMS of the median charge)
    std::vector<size_t> GetMIPindices(const std::vector<double>& dQ,
				      const double& median,
				      const double& rms,
				      const double& alpha);

    // get a subvector containing the entries for the list of indices provided
    std::vector<double> GetSubVector(const std::vector<double> v,
				     const std::vector<size_t> idx);
    
    // get the median value of a vector
    double GetMedian(std::vector<double>& v);

    // get average of values
    double GetAvg(const std::vector<double>& v);
    // get rms of values
    double GetRms(const std::vector<double>& v, const double avg);
    double GetRms(const std::vector<double>& v);
    double GetRms(const std::vector<double>& MIPdS,
		  const std::vector<double>& MIPdQ,
		  const double& slope,
		  const double& intercept);

    // get max index and value
    std::pair<size_t,double> GetMaxIndex(const std::vector<double>& v);

    // function to calculate linear fit to list of points
    std::pair<double,double> GetLinearFit(const std::vector<double>& x,
					  const std::vector<double>& y);


    // function to calculate Bragg Peak area under the curve (taking into account best fit to MIP slope)
    double GetBraggArea(const std::vector<double>& dS, const std::vector<double>& dQ,
			const size_t& MIPendIdx, const size_t& braggIdx,
			const double& MIPm, const double& MIPs);

    // tree
    TTree* _tree;
    // final muon energy
    double _Estop;
    // truncated charge vector
    std::vector<double> _dQ;
    // dS vector
    std::vector<double> _dS;    
    // subset of hits in MIP region
    std::vector<double> _MIPdS;
    std::vector<double> _MIPdQ;
    std::vector<double> _muonMIPdS;
    std::vector<double> _muonMIPdQ;
    
    // intercept obtained by fitting MIP points to a line
    double _MIPm;
    // slope obtained by fitting MIP points to a line
    double _MIPs;
    // median hit charge value of all MIP points 15 cm before the bragg peak
    double _MIPmedian;
    // rms hit charge value from all MIP points 15 cm before the bragg peak
    double _MIPrms;
    // rms as above, but calculated after correcting for the slope, and only for those hits within a certain RMS range of the median
    double _MIPrms_corr;
    // index of bragg peak
    int _braggIdx;
    // bragg peak charge
    double _braggQ;
    // expected bragg amp (i.e. continuing MIP fit 'till braggIndx)
    double _braggExpected;
    // bragg amplitude
    double _braggAmp;
    // bragg peak area
    double _braggArea;
    // distance in cm between bragg peak and end of track
    double _braggDistToEnd;
    // MIP end position (15 cm before Bragg)
    int _MIPendIdx;
    

    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
