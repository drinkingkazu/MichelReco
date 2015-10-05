import sys, os

import ROOT
from ROOT import larlite

import matplotlib.pyplot as plt
plt.ion()
import numpy as np

from ROOT import larutil
t2cm = larutil.GeometryHelper.GetME().TimeToCm()
w2cm = larutil.GeometryHelper.GetME().WireToCm()

manager = larlite.storage_manager()

manager.reset()
for i in xrange(len(sys.argv)-1):
    manager.add_in_filename(sys.argv[i+1])
manager.set_io_mode(larlite.storage_manager.kREAD)
manager.open()

fig, ax = plt.subplots(nrows=1,figsize=(10,7))

colors = ['b','c','m','g','k','y']

while manager.next_event():

    #muons   = manager.get_data(larlite.data.kCluster,'rawcluster')
    muons   = manager.get_data(larlite.data.kCluster,'muon')
    michels = manager.get_data(larlite.data.kCluster,'michel')
    opflash = manager.get_data(larlite.data.kOpFlash,'opflash')
    hits    = manager.get_data(larlite.data.kHit,'gaushit')
    muon_ass = manager.get_data(larlite.data.kAssociation,'muon')
    michel_ass = manager.get_data(larlite.data.kAssociation,'michel')
    muon_index_v = muon_ass.association(muons.id(), hits.id())
    michel_index_v = michel_ass.association(michels.id(), hits.id())

    if (muon_index_v.size() == 0):
        continue

    # save michel cluster
    michel_z = []
    michel_x = []
    

    # we have a cluster -> get the average z-position and total charge
    all_clus_z_avg = []
    all_clus_q_tot = []
    all_clus_z = []
    all_clus_x = []

    for clus in xrange(len(muon_index_v)):

        # skip planes 0 and 1
        if (hits[muon_index_v[clus][0]].WireID().Plane != 2):
            continue

        clus_z_avg = 0.
        clus_q_tot = 0.
        
        clus_z = []
        clus_x = []
        min_x = 1036.
        
        for hit in muon_index_v[clus]:
            hit_z = hits[hit].WireID().Wire*w2cm
            hit_x = hits[hit].PeakTime()*t2cm
            if (hit_x < min_x):
                min_x = hit_x
            clus_z.append(hit_z)
            clus_x.append(hit_x)
            clus_z_avg += hit_z
            clus_q_tot += hits[hit].Integral()

        clus_z_avg /= muon_index_v[clus].size()

        clus_z = np.array(clus_z)
        clus_x = np.array(clus_x)-min_x+10

        all_clus_z.append(clus_z)
        all_clus_x.append(clus_x)
        all_clus_z_avg.append(clus_z_avg)
        all_clus_q_tot.append(clus_q_tot)


    # print the number of clusters we have
    print 'numbr of clusters : %i'%len(all_clus_z)

    # now get all optical flashes
    print 'number of flashes : %i'%opflash.size()

    # keep track of all flash z positions and their charge
    flash_z  = []
    flash_pe = []
    for n in xrange(opflash.size()):
        
        flash = opflash.at(n)

        flash_z.append(flash.ZCenter())
        flash_pe.append(flash.TotalPE())

    flash_y = np.ones(len(flash_z))

    ax.clear()
    fig.gca()

    plt.scatter(flash_z,flash_y,s=flash_pe,c='r')
    for i in xrange(len(all_clus_z)):
        plt.axvline(all_clus_z_avg[i],lw=2,color=colors[i%6])
        plt.scatter(all_clus_z_avg[i],[256.],s=200,c=colors[i%6])
        plt.scatter(all_clus_z[i],all_clus_x[i],s=20,c=colors[i%6],edgecolor='None')
    plt.grid()
    plt.show()

    usrinput = raw_input("Hit Enter: next evt  ||  q: exit viewer\n")
    if ( usrinput == "q" ):
        sys.exit(0)
