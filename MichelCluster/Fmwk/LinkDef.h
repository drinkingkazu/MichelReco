//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace michel+;
#pragma link C++ namespace michel::msg+;
#pragma link C++ class michel::HitPt+;
#pragma link C++ class michel::Michel+;
#pragma link C++ class michel::MichelCluster+;
#pragma link C++ class std::vector<michel::MichelCluster>+;
#pragma link C++ class michel::MichelException+;
#pragma link C++ class michel::MichelReco+;
#pragma link C++ class michel::MichelAnaBase+;
#pragma link C++ class michel::BaseMichelAlgo+;
#pragma link C++ class michel::BaseAlgMerger+;
#pragma link C++ class michel::BaseAlgBinaryMerger+;
#pragma link C++ class michel::BaseAlgIdentifier+;
#pragma link C++ class michel::BaseAlgBoundary+;

#pragma link C++ class cmtool::CMergeBookKeeper+;
#pragma link C++ class cmtool::CMTException+;

#pragma link C++ class larlite::MichelRecoDriver+;
//ADD_NEW_CLASS ... do not change this line
#endif













