import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from ROOT import gSystem,TMath
from ROOT import michel
from ROOT import larlite as fmwk
from ROOT import larutil


# Create ana_processor instance
my_proc = fmwk.ana_processor()

my_proc.enable_filter(True)
the_filter = fmwk.StoppingFilter()
the_filter.Filter(True)
the_filter.Stopping(True)
the_filter.Pure(True)

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify analysis output root file name
my_proc.set_ana_output_file("out.root");

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

stopping = fmwk.TagStoppingMuon()
stopping.SetClusterProducer('cccluster')
stopping.SetAlgo(michel.TruncatedQBoundary())

#my_proc.add_process(the_filter)
my_proc.add_process(stopping)

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit(0)

