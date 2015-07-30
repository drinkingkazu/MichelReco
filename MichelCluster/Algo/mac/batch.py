import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

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
my_unit.SetClusterProducer("fuzzycluster")

# Get manager for michel reco
mgr = my_unit.Algo()

# Attach algorithm for merging
mgr.SetAlgo(michel.kClusterMerger, michel.ToyMerger())

# Attach algorithm for boundary finding
mgr.SetAlgo(michel.kBoundaryFinder, michel.ToyBoundary())

# Attach algorithm for finding michel cluster
mgr.SetAlgo(michel.kMichelID, michel.ToyMichelID())

# Attach algorithm to recluster michel
mgr.SetAlgo(michel.kMichelCluster, michel.ToyMichelCluster())

# Attach ana unit
my_proc.add_process(my_unit)

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run();

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
