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

my_proc.set_ana_output_file("ana.root")

# Michel reco driver code
my_unit = fmwk.MichelRecoDriver()
#my_unit.SetClusterProducer("fuzzycluster")
my_unit.SetClusterProducer("linecluster")

# Get manager for michel reco
mgr = my_unit.GetManager()

#mgr.SetVerbosity(michel.msg.kDEBUG)

# Attach algorithm for merging
mgr.SetAlgo(michel.kClusterMerger, 
            michel.EdgeMerger())

# Attach algorithm for boundary finding
mgr.SetAlgo(michel.kBoundaryFinder, 
            michel.TruncatedQBoundary())

# Attach algorithm for finding michel cluster
mgr.SetAlgo(michel.kMichelID, 
            michel.ForwardMichelID())

# Attach algorithm to recluster michel
mgr.SetAlgo(michel.kMichelCluster, 
            michel.RadiusMichelCluster())

# Attach ana unit

mgr.AddAna(michel.CosmicAna())

# add process to get moving
my_proc.add_process(my_unit)


#Write aho unit to get out MC vars whatever

# aho_ana=fmwk.AhoAna()
# my_proc.add_process(aho_ana)

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
#my_proc.run(145,5);

while (my_proc.process_event()):

    print "next event..."
    alg = my_unit.Algo()
    michels = alg.GetResult()
    for x in xrange(michels.size()):
        m = michels[x]._michel
        print "size of this michel: ",  m.size()
    print

#my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
