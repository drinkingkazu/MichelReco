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

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
#my_proc.set_io_mode(fmwk.storage_manager.kREAD)

#my_proc.set_ana_output_file("michel_tree.root")

my_proc.set_output_file("michel_clusters.root")

#Cheat and send in signal
my_proc.enable_filter(True)
the_filter = fmwk.MichelFilter()
#the_filter = fmwk.RemoveMichel()

# Michel reco driver code
my_unit = fmwk.MichelRecoDriver()
my_unit.SetClusterProducer("fuzzycluster")
my_unit.SetUseMC(True)
#my_unit.SetClusterProducer("rawcluster")
my_unit.SetEField(0.5)

# set here if you want to save michels as an output cluster
my_unit.saveMichelClusters(True)

# Get manager for michel reco
mgr = my_unit.GetManager()

#mgr.SetVerbosity(michel.msg.kDEBUG)
#mgr.SetDebug(True)

# Attach algorithm for merging
mgr.AddMergingAlgo(michel.EdgeMerger())

# Attach algorithm for boundary finding
truncBound = michel.TruncatedQBoundary()
truncBound.SetMaxDistanceTruncatedPeaks(5)
chiBound   = michel.ChiBoundary()
matchBound = michel.MatchBoundaries()
matchBound.SetMaxDistanceTruncatedPeaks(10)
matchBound.SetMaxCovarianceAtStart(0.8)
covariance = michel.CovarianceFollowBoundary()
covariance.SetMaxDistanceTruncatedPeaks(10)
mgr.AddAlgo(covariance)

# Attach algorithm for finding michel cluster
findMichel = michel.ForwardMichelID()
findMichel.SetMaxMichelHits(0)
mgr.AddAlgo(findMichel)

# MID finding algorithm
midalgo = michel.DecideIfStoppingMuon()
midalgo.SetChiMin(0.9)
midalgo.SetFracMinHits(0.7)
midalgo.SetHitRadius(30)
midalgo.SetMaxDist(3.0)
midalgo.SetMinBadHits(10)
#midalgo.SetVerbosity(michel.msg.kDEBUG)
mgr.AddAlgo(midalgo)
minlength = michel.CutOnMuonLength()
minlength.SetMinMuonLength(10)
mgr.AddAlgo(minlength)
minlinearity = michel.CutOnMuonLinearity()
minlinearity.SetChiMin(0.8)
minlinearity.SetFracMinHits(0.5)
mgr.AddAlgo(minlinearity)

# Attach algorithm to recluster michel
supersonic = michel.SuperSonicClusterer()
stepsonic  = michel.StepSuperSonicCluster()
mgr.AddAlgo(stepsonic)

# Attach ana unit
mgr.AddAna(michel.CosmicAna())

# add process to get moving
#my_proc.add_process(the_filter)
my_proc.add_process(my_unit)

my_proc.set_data_to_write(fmwk.data.kHit,'cchit')
#my_proc.set_data_to_write(fmwk.data.kHit,'gaushit')
my_proc.set_data_to_write(fmwk.data.kCluster,'michel')
#my_proc.set_data_to_write(fmwk.data.kCluster,'rawclusters')
my_proc.set_data_to_write(fmwk.data.kAssociation,'michel')
#my_proc.set_data_to_write(fmwk.data.kAssociation,'rawclusters')
#Write aho unit to get out MC vars whatever

#my_proc.enable_event_alignment(False)

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
#my_proc.run(145,5);
my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
