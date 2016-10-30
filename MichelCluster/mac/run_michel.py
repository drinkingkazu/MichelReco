import sys
import ROOT
from larlite import larlite as fmwk
from ROOT import michel
from highStats_algo_devel import PrepareMichelAlgo
import time
import os

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

my_proc.set_ana_output_file(  "michel_tree.root"  )
my_proc.set_output_file    ("michel_clusters.root")

my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

#########################
# Michel reco driver code
my_unit = fmwk.MichelRecoDriver()
#my_unit.SetClusterProducer("fuzzycluster")
my_unit.SetClusterProducer("pandoraCosmic")
#my_unit.SetClusterProducer("linecluster")
my_unit.saveMichelClusters(True)
#my_unit.SetUseMC(False)
my_unit.SetMinClusSize(15)
my_unit.SetAvgGain(5.3)
my_unit.SetNChannels(8256)
my_unit.set_outtxt_file( "michel_eventinfo.txt" )

# load ch-by-ch gain values
chchgains_f = open('%s/MichelReco/MichelCluster/mac/ChChGain.txt'%(os.environ['LARLITE_USERDEVDIR']),'r')
for line in chchgains_f:
    words = line.split()
    chan = int(words[0])
    gain = float(words[1])
    my_unit.SetChGain(chan,5.3)

###########################################################
# set here if you want to save michels as an output cluster
#my_unit.saveMichelClusters(True)

#############################
# Get manager for michel reco
mgr = my_unit.GetManager()

#####################################
# Debug options
# mgr.SetVerbosity(michel.msg.kDEBUG)
# mgr.SetDebug(True)

#############################
# Attach cluster filter algo
filterAlgo =  michel.FilterStraightLineClusters()
filterAlgo.setMinRMS(0.5)
mgr.AddFilteringAlgo(filterAlgo)

##############################
# Attach algorithm for merging
mgr.AddMergingAlgo(michel.EdgeMerger())


algoList = PrepareMichelAlgo()

for algo in algoList:
    mgr.AddAlgo(algo)

# Attach ana unit
mgr.AddAna(michel.CosmicAna())

# add process

my_unit.FilterEvents(False)

my_proc.add_process(my_unit)

filterMichels = fmwk.FilterMichels()
filterMichels.SetProducer("electron")
my_proc.add_process(filterMichels)

reco = fmwk.RecoEffStudy()
reco.SetDebug(False)
reco.SetMCShowerProducer("mcreco2")
reco.SetGoodMatchDistance(30)
my_proc.add_process(reco)

my_proc.enable_filter(True)

#my_proc.set_data_to_write(fmwk.data.kHit,'cchit')
#my_proc.set_data_to_write(fmwk.data.kHit,'linecluster')
my_proc.set_data_to_write(fmwk.data.kHit,'gaushit')
my_proc.set_data_to_write(fmwk.data.kCluster,'electron')
my_proc.set_data_to_write(fmwk.data.kCluster,'photon')
my_proc.set_data_to_write(fmwk.data.kAssociation,'electron')
my_proc.set_data_to_write(fmwk.data.kAssociation,'photon')


print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
#my_proc.run(0,100);
now = time.clock()
my_proc.run()
dt = time.clock() - now

print 'Total running time : %i [s]'%(int(dt))

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
