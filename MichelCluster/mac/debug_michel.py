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


algos_do = []
algos_test = []

#########################################
# calculate various cluster parameters...
ctrunk = michel.CalcTruncated()
#ctrunk.SetCovarianceWindowSize(int s)      
#ctrunk.SetTruncatedQWindowSize(int s)      
#ctrunk.SetPAbove(double p)                 
#ctrunk.SetMinWindowSize(int w)             
#ctrunk.SetEdgeEffectFix(int e)             
mgr.AddAlgo(ctrunk)
algos_do.append(ctrunk)

########################################
# Attach algorithm for boundary finding
braggalgo = michel.FindBraggPeak()
braggalgo.SetMinBraggArea(200)
#mgr.AddAlgo(braggalgo)
#algos_test.append(braggalgo)
boundaryalgo = michel.BoundaryFromTQMaxQ()
boundaryalgo.SetMaxDistancesTruncatedQMaxQ(15)
boundaryalgo.SetVerbosity(michel.msg.kDEBUG)
algos_test.append(boundaryalgo)

##############################################
# Attach algo for charge spectrum requirements
closepeaks = michel.RequireCloseTruncatedPeaks()
closepeaks.SetMaxDistanceTruncatedPeaks(5)
#mgr.AddAlgo(closepeaks)

covdip = michel.RequireCovarianceDip()
covdip.SetCovarianceDipCutoff(0.9)
#mgr.AddAlgo(covdip)

slopeflip = michel.RequireSlopeSignFlip()
#mgr.AddAlgo(slopeflip)

lowcovbound = michel.RequireBoundaryInLowCov()
lowcovbound.SetMaxCovarianceAtStart(0.8)
lowcovbound.SetVerbosity(michel.msg.kDEBUG)
algos_test.append(lowcovbound)

#############################################
# Attach algorithm for finding michel cluster
findMichel = michel.ForwardMichelID()
findMichel.SetMaxMichelHits(0)
findMichel.SetVerbosity(michel.msg.kDEBUG)
algos_test.append(findMichel)

#########################
# MID finding algorithms
midalgo = michel.DecideIfStoppingMuon()
midalgo.SetChiMin     ( 0.9 )
midalgo.SetFracMinHits( 0.7 )
midalgo.SetHitRadius  ( 30  )
midalgo.SetMaxDist    ( 3.0 )
midalgo.SetMinBadHits ( 10  )
midalgo.SetVerbosity(michel.msg.kDEBUG)
algos_test.append(midalgo)

minlength = michel.CutOnMuonLength()
minlength.SetMinMuonLength(10)
algos_test.append(minlength)

minlinearity = michel.CutOnMuonLinearity()
minlinearity.SetChiMin     ( 0.8 )
minlinearity.SetFracMinHits( 0.5 )
algos_test.append(minlinearity)

#########################################################
# MID filter that removes michels close to wire gaps/edges
fidvolfilter = michel.CutOnFiducialVolume()
import parse_fiducial_volume_definitions as fidparser

wires_to_exclude_min, \
    wires_to_exclude_max, \
    times_to_exclude_min, \
    times_to_exclude_max = fidparser.list_wires_times_to_exclude()

fidvolfilter.SetExcludedWireRanges(wires_to_exclude_min,wires_to_exclude_max)
fidvolfilter.SetExcludedTimeRanges(times_to_exclude_min,times_to_exclude_max)
algos_test.append(fidvolfilter)

######################################
# Attach algorithm to recluster michel
supersonic = michel.SuperSonicClusterer()
supersonic.SetMergeTillConverge(True)
supersonic.SetUseHitRadius     (True)
supersonic.SetMaxRadius( 15 )
supersonic.SetHitRadius(  3 )
algos_test.append(supersonic)

#########################
# cone-finding algorithm
conefinder = michel.ConeHitFinder()
conefinder.SetMaxRadius               ( 20 )
conefinder.SetMaxPerpendicularDistance(  3 )
algos_test.append(conefinder)

##############################################
# final mid algo cutting on num of michel hits
michelhits = michel.CutOnMichelNumHits()
michelhits.SetMinMichelHits (  5 )
michelhits.SetMaxMichelHits ( 35 )
algos_test.append(michelhits)

#########################################
# remove weird horizontal tracks from PMT
pmtremoved = michel.RemoveFakePMTSignals()
pmtremoved.SetMaxErrorTime(0.1)
algos_test.append(pmtremoved)

#############################################
# require large angle between michel and muon
largeangle = michel.RequireLargeAngle()
largeangle.SetMinAngle(30.*3.14/180.)
largeangle.SetMinStraightMichelHits(5)
algos_test.append(largeangle)

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

        print
        print 'processing cluster %i'%c_idx
        proceed = True
        for algo in algos_do:
            res = algo.ProcessCluster(clus,hit_v)
            if (res == False):
                print '-> Stopping @ Algo %s'%algo.Name()
                proceed = False
                break
        if (proceed == False):
            continue

        dS = np.array(clus._s_v)
        dQ = np.array(clus._t_mean_v)
        dX = np.abs(np.array(clus._chi2_v))
        if (np.amax(dQ) < 50):
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
        axarr[1].plot(dS,hit_q,'o--',color='b')
        axarr[1].plot(dS,dQ,'o--',color='c')
        if (clus._boundary > 0 and clus._boundary < len(dS)):
            axarr[1].scatter(dS[clus._boundary],dQ[clus._boundary],c='k',s=50)
        axarr[1].grid()
        ax_chi.plot(dS,dX,'o--',color='r')
        for ti in ax_chi.get_yticklabels():
            ti.set_color('r')

        for algo in algos_test:
            res = algo.ProcessCluster(clus,hit_v)
            if (res == False):
                print '-> Stopping @ Algo %s'%algo.Name()
                proceed = False
                break

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