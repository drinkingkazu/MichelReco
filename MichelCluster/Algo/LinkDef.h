//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;



/// ClusterMerger
#pragma link C++ class michel::EdgeMerger+;

/// BoundaryFinder
#pragma link C++ class michel::ChiBoundary+;
#pragma link C++ class michel::TruncatedQBoundary+;
#pragma link C++ class michel::TSpectrumBoundary+;
#pragma link C++ class michel::MatchBoundaries+;
#pragma link C++ class michel::CovarianceFollowBoundary+;

/// MichelID
#pragma link C++ class michel::ForwardMichelID+;

/// MichelCluster
#pragma link C++ class michel::RadiusMichelCluster+;
#pragma link C++ class michel::SuperSonicClusterer+;
#pragma link C++ class michel::StepAroundCluster+;
#pragma link C++ class michel::StepSuperSonicCluster+;

/// MIDFilter
#pragma link C++ class michel::DecideIfStoppingMuon+;
#pragma link C++ class michel::CutOnMuonLength+;
#pragma link C++ class michel::CutOnMuonLinearity+;

/// Various method classes
#pragma link C++ class michel::ClusterVectorCalculator+;

// #pragma link C++ class michel::ToyMerger+;
// #pragma link C++ class michel::ToyBoundary+;
// #pragma link C++ class michel::ToyMichelID+;
// #pragma link C++ class michel::ToyMichelCluster+;

//ADD_NEW_CLASS ... do not change this line
#endif
