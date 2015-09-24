import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

import ROOT
ROOT.gSystem.Load("libMichelReco_MichelClusterFmwk")
from larlite import larlite as fmwk
from ROOT import michel
import matplotlib.pyplot as plt
import numpy as np
from highStats_algo import PrepareMichelAlgo

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

my_proc.set_io_mode(fmwk.storage_manager.kREAD)

my_proc.set_ana_output_file(  "michel_tree.root"  )
my_proc.set_output_file    ("michel_clusters.root")

#########################
# Michel reco driver code
my_unit = fmwk.MichelRecoDriver()
#my_unit.SetClusterProducer("fuzzycluster")
#my_unit.SetClusterProducer("linecluster")
my_unit.SetClusterProducer("rawcluster")
my_unit.SetUseMC(False)
my_unit.SetEField(0.27)
my_unit.SetMinClusSize(15)

###########################################################
# set here if you want to save michels as an output cluster
# my_unit.saveMichelClusters(True)

#############################
# Get manager for michel reco
mgr = my_unit.GetManager()

#####################################
# Debug options
# mgr.SetVerbosity(michel.msg.kDEBUG)
# mgr.SetDebug(True)

##############################
# Attach algorithm for merging
mgr.AddMergingAlgo(michel.EdgeMerger())


algoList = PrepareMichelAlgo()
for algo in algoList:
    mgr.AddAlgo(algo)

# add process
my_proc.add_process(my_unit)

my_proc.set_data_to_write(fmwk.data.kHit,'cchit')
#my_proc.set_data_to_write(fmwk.data.kHit,'gaushit')
my_proc.set_data_to_write(fmwk.data.kCluster,'michel')
#my_proc.set_data_to_write(fmwk.data.kCluster,'rawclusters')
my_proc.set_data_to_write(fmwk.data.kAssociation,'michel')
#my_proc.set_data_to_write(fmwk.data.kAssociation,'rawclusters')


print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
#my_proc.run(0,100);

# get algorithm chain
algo_chain = mgr.GetAlgoVector()

plt.ion()
fig, axarr = plt.subplots(2,1,figsize=(12,8))
ax_chi = axarr[1].twinx()

while my_proc.process_event():

    hit_v   = mgr.GetHitVector()
    all_hit_w = []
    all_hit_t = []
    for h in hit_v:
        if (h._pl == 2):
            all_hit_w.append(h._w)
            all_hit_t.append(h._t)

    in_clus = mgr.GetMergedClusters()
    print 'There are %i'%len(in_clus)
    # for each cluster, apply the full algo chain
    for c_idx,clus in enumerate(in_clus):

        axarr[0].cla()
        axarr[1].cla()
        ax_chi.cla()
        fig.gca()

        proceed = True

        for algo in algoList:
            res = algo.ProcessCluster(clus,hit_v)
            if (res == False):
                #print '-> Stopping @ Algo %s'%algo.Name()
                proceed = False
                break

        dS = np.array(clus._s_v)
        dQ = np.array(clus._t_mean_v)
        dX = np.abs(np.array(clus._chi2_v))
        if (np.amax(dQ) < 50):
            continue


        if (proceed == False):
            continue

        hits = clus._hits;
        hit_t = []
        hit_w = []
        hit_q = []
        a_hits = clus._all_hits;
        a_hit_t = []
        a_hit_w = []
        a_hit_q = []
        for h in hits:
            hit_t.append(h._t)
            hit_w.append(h._w)
            hit_q.append(h._q)
        hit_t = np.array(hit_t)
        hit_w = np.array(hit_w)
        hit_q = np.array(hit_q)
        for h in a_hits:
            a_hit_t.append(h._t)
            a_hit_w.append(h._w)
            a_hit_q.append(h._q)
        a_hit_t = np.array(a_hit_t)
        a_hit_w = np.array(a_hit_w)
        a_hit_q = np.array(a_hit_q)

        # get michel cluster w,t range
        w_max = np.amax(hit_w)
        w_min = np.amin(hit_w)
        t_max = np.amax(hit_t)
        t_min = np.amin(hit_t)

        axarr[0].scatter(all_hit_w,all_hit_t,s=60,edgecolor='k',facecolor='none')
        axarr[0].scatter(a_hit_w,a_hit_t,s=60,edgecolor='none',facecolor='r')
        axarr[0].scatter(hit_w,hit_t,c=hit_q,s=60,edgecolor='none')

        axarr[0].set_xlim([w_min-20,w_max+20])
        axarr[0].set_ylim([t_min-20,t_max+20])
        axarr[0].grid()
        bound = clus.GetBoundary()
        if (bound > 0 and bound < len(dS)):
            axarr[1].axvline(dS[bound],c='k',lw=3)
        axarr[1].plot(dS,hit_q,'o--',color='b')
        axarr[1].plot(dS,dQ,'o--',color='c')
        axarr[1].grid()
        axarr[1].set_ylim([0,np.amax(dQ)])
        ax_chi.plot(dS,dX,'o--',color='r')
        for ti in ax_chi.get_yticklabels():
            ti.set_color('r')


        # michel hits
        michel_t = []
        michel_w = []
        for h in clus._michel:
            michel_t.append(h._t)
            michel_w.append(h._w)
        michel_t = np.array(michel_t)
        michel_w = np.array(michel_w)
        if (len(michel_w) > 0):
            axarr[0].scatter(michel_w,michel_t,c='k',s=60,edgecolor='none')

        fig.canvas
        fig.canvas.draw()

        usrinput = raw_input("Hit Enter: next evt  || int: go to event number ||  q: exit viewer\n")
        if ( usrinput == "q" ):
            sys.exit(0)


        

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
