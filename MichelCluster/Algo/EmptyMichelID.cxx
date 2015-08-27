#ifndef MICHELCLUSTER_EMPTYMICHELID_CXX
#define MICHELCLUSTER_EMPTYMICHELID_CXX

#include "EmptyMichelID.h"
#include "Fmwk/MichelException.h"

namespace michel {

Michel EmptyMichelID::Identify(const MichelCluster& cluster, bool& forward)
{

	Michel electron;

	for ( size_t i = 0 ; i < cluster._ordered_pts.size(); ++i) {
		electron.push_back(cluster._hits[cluster._ordered_pts[i]]);
	}


	return electron;
}

}
#endif
