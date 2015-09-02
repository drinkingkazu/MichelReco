#ifndef MICHELCLUSTER_EMPTYMICHELID_CXX
#define MICHELCLUSTER_EMPTYMICHELID_CXX

#include "EmptyMichelID.h"
#include "Fmwk/MichelException.h"

namespace michel {

void EmptyMichelID::ProcessCluster(MichelCluster& cluster,
				   const std::vector<HitPt>& hits)
{
	auto electron = cluster._michel;
	for ( size_t i = 0 ; i < cluster._hits.size(); ++i) 
		electron.push_back(cluster._hits[i]);
	return;
}

}
#endif
