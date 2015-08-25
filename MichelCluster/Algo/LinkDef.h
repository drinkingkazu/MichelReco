//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class michel::ToyMerger+;
#pragma link C++ class michel::ToyBoundary+;
#pragma link C++ class michel::ChiBoundary+;
#pragma link C++ class michel::TruncatedQBoundary+;
#pragma link C++ class michel::TSpectrumBoundary+;
#pragma link C++ class michel::ToyMichelID+;
#pragma link C++ class michel::ToyMichelCluster+;
#pragma link C++ class michel::EdgeMerger+;
#pragma link C++ class michel::ForwardMichelID+;
#pragma link C++ class michel::RadiusMichelCluster+;
#pragma link C++ class michel::SuperSonicClusterer+;
#pragma link C++ class michel::ClusterVectorCalculator+;

//ADD_NEW_CLASS ... do not change this line
#endif
