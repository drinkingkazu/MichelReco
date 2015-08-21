#ifndef LARLITE_MICHELCLUSTERUNITTEST_CXX
#define LARLITE_MICHELCLUSTERUNITTEST_CXX

#include "MichelClusterUnitTest.h"

namespace larlite {

bool MichelClusterUnitTest::initialize() {

    return true;
}

bool MichelClusterUnitTest::analyze(storage_manager* storage) {

    //Grab the clusters
    auto ev_cluster = storage->get_data<event_cluster>("michel");
    if (!ev_cluster) {
        print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, cluster!"));
        return false;
    }

    if (!ev_cluster->size())
        return false;

    event_hit* ev_hit = nullptr;
    auto const& ass_hit_v = storage->find_one_ass(ev_cluster->id(), ev_hit, ev_cluster->name());

    if (!ev_hit) {
        print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find hits."));
        return false;
    }
    if (!ass_hit_v.size())
        print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find any hit associations."));

    //Loop over clusters
    for (size_t i = 0; i < ev_cluster->size(); ++i) {

        if (!ass_hit_v[i].size())
            continue;

        for (auto const& index : ass_hit_v[i])
            std::cout << "ass_hit_v[" << i << "] tells me to look in ev_hit at index " << index << std::endl;

        std::cout << "Finished loop over over associated hits for this cluster!" << std::endl;

    }
    return true;
}

bool MichelClusterUnitTest::finalize() {
    return true;
}

}
#endif
