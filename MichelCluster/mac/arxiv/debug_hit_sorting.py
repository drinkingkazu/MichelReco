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
#mgr.AddMergingAlgo(michel.EdgeMerger())


#########################################
# calculate various cluster parameters...
ctrunk = michel.CalcTruncated()
mgr.AddAlgo(ctrunk)


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
fig, ax =  plt.subplots(figsize=(12,8))

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

        ax.cla()
        fig.gca()

        print
        print 'processing cluster %i'%c_idx


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
        w_max = np.amax(a_hit_w)
        w_min = np.amin(a_hit_w)
        t_max = np.amax(a_hit_t)
        t_min = np.amin(a_hit_t)

        ax.scatter(all_hit_w,all_hit_t,s=60,edgecolor='k',facecolor='none')
        ax.scatter(a_hit_w,a_hit_t,s=60,edgecolor='r',facecolor='r')
        ax.scatter(hit_w,hit_t,c=hit_q,s=60,edgecolor='none')
        ax.set_xlim([w_min-10,w_max+10])
        ax.set_ylim([t_min-10,t_max+10])
        ax.grid()

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
