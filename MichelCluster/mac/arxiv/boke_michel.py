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
my_unit.SetEField(0.5)

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

#########################################
# calculate various cluster parameters...
ctrunk = michel.CalcTruncated()
#ctrunk.SetCovarianceWindowSize(int s)      
#ctrunk.SetTruncatedQWindowSize(int s)      
#ctrunk.SetPAbove(double p)                 
#ctrunk.SetMinWindowSize(int w)             
#ctrunk.SetEdgeEffectFix(int e)             
mgr.AddAlgo(ctrunk)

########################################
# Attach algorithm for boundary finding
boundaryalgo = michel.BoundaryFromTQMaxQ()
boundaryalgo.SetMaxDistancesTruncatedQMaxQ(5)
mgr.AddAlgo(boundaryalgo)

##############################################
# Attach algo for charge spectrum requirements
closepeaks = michel.RequireCloseTruncatedPeaks()
closepeaks.SetMaxDistanceTruncatedPeaks(5)
mgr.AddAlgo(closepeaks)

covdip = michel.RequireCovarianceDip()
covdip.SetCovarianceDipCutoff(0.9)
mgr.AddAlgo(covdip)

slopeflip = michel.RequireSlopeSignFlip()
mgr.AddAlgo(slopeflip)

lowcovbound = michel.RequireBoundaryInLowCov()
lowcovbound.SetMaxCovarianceAtStart(0.8)
mgr.AddAlgo(lowcovbound)

#############################################
# Attach algorithm for finding michel cluster
findMichel = michel.ForwardMichelID()
findMichel.SetMaxMichelHits(0)
mgr.AddAlgo(findMichel)

#########################
# MID finding algorithms
midalgo = michel.DecideIfStoppingMuon()
midalgo.SetChiMin     ( 0.9 )
midalgo.SetFracMinHits( 0.7 )
midalgo.SetHitRadius  ( 30  )
midalgo.SetMaxDist    ( 3.0 )
midalgo.SetMinBadHits ( 10  )
mgr.AddAlgo(midalgo)

minlength = michel.CutOnMuonLength()
minlength.SetMinMuonLength(10)
mgr.AddAlgo(minlength)

minlinearity = michel.CutOnMuonLinearity()
minlinearity.SetChiMin     ( 0.8 )
minlinearity.SetFracMinHits( 0.5 )
mgr.AddAlgo(minlinearity)

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
mgr.AddAlgo(fidvolfilter)

######################################
# Attach algorithm to recluster michel
supersonic = michel.SuperSonicClusterer()
supersonic.SetMergeTillConverge(True)
supersonic.SetUseHitRadius     (True)
supersonic.SetMaxRadius( 15 )
supersonic.SetHitRadius(  3 )
mgr.AddAlgo(supersonic)

#########################
# cone-finding algorithm
conefinder = michel.ConeHitFinder()
conefinder.SetMaxRadius               ( 20 )
conefinder.SetMaxPerpendicularDistance(  3 )
#mgr.AddAlgo(conefinder)

##############################################
# final mid algo cutting on num of michel hits
michelhits = michel.CutOnMichelNumHits()
michelhits.SetMinMichelHits (  5 )
michelhits.SetMaxMichelHits ( 35 )
mgr.AddAlgo(michelhits)

#########################################
# remove weird horizontal tracks from PMT
pmtremoved = michel.RemoveFakePMTSignals()
pmtremoved.SetMaxErrorTime(0.1)
mgr.AddAlgo(pmtremoved)

#############################################
# require large angle between michel and muon
largeangle = michel.RequireLargeAngle()
largeangle.SetMinAngle(30.*3.14/180.)
largeangle.SetMinStraightMichelHits(5)
mgr.AddAlgo(largeangle)

# Attach ana unit
mgr.AddAna(michel.CosmicAna())

# add process
my_proc.add_process(my_unit)

#my_proc.set_data_to_write(fmwk.data.kHit,'cchit')
my_proc.set_data_to_write(fmwk.data.kHit,'gaushit')
my_proc.set_data_to_write(fmwk.data.kCluster,'michel')
#my_proc.set_data_to_write(fmwk.data.kCluster,'rawclusters')
my_proc.set_data_to_write(fmwk.data.kAssociation,'michel')
#my_proc.set_data_to_write(fmwk.data.kAssociation,'rawclusters')


print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
#my_proc.run(0,100);
my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
