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
from highStats_algo import PrepareMichelAlgo

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

my_proc.set_ana_output_file(  "michel_tree.root"  )
my_proc.set_output_file    ("michel_clusters.root")

#########################
# Michel reco driver code
my_unit = fmwk.MichelRecoDriver()
#my_unit.SetClusterProducer("fuzzycluster")
#my_unit.SetClusterProducer("rawcluster")
my_unit.SetClusterProducer("pandoraCosmic")
my_unit.saveMichelClusters(True)
#my_unit.SetClusterProducer("linecluster")
my_unit.SetUseMC(False)
my_unit.SetEField(0.5)
my_unit.SetMinClusSize(15)

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
my_proc.add_process(my_unit)

#my_proc.set_data_to_write(fmwk.data.kHit,'cchit')
#my_proc.set_data_to_write(fmwk.data.kHit,'linecluster')
#my_proc.set_data_to_write(fmwk.data.kHit,'gaushit')
my_proc.set_data_to_write(fmwk.data.kCluster,'michel')
my_proc.set_data_to_write(fmwk.data.kCluster,'muon')
#my_proc.set_data_to_write(fmwk.data.kCluster,'rawclusters')
my_proc.set_data_to_write(fmwk.data.kAssociation,'michel')
my_proc.set_data_to_write(fmwk.data.kAssociation,'muon')
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
