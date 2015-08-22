# Load libraries
import ROOT, sys, os
from ROOT import *
import numpy as np
gSystem.Load('libRecoTool_CMToolApp.so')
# Now import ana_processor & your class. For this example, ana_base.
if len(sys.argv) < 2:
    print
    print "*** Improper usage. Usage: python %s /path/to/input/file.root ***" % sys.argv[0]
    print

handscan_textfilename = os.environ['LARLITE_USERDEVDIR']+'/MichelReco/MichelCluster/mac/michel_handscanning_results.txt'

#dict of (run, subrun, eventid) => [list of cluster indexes scanned]
#this dictionary keeps track of what you have already scanned
scanned_events_dict = {}

my_proc = larlite.ana_processor()
my_proc.set_verbosity(larlite.msg.kDEBUG)

my_proc.set_io_mode(larlite.storage_manager.kREAD)

#args should be input file name
for x in xrange(len(sys.argv)-1):

    my_proc.add_input_file(sys.argv[x+1])

my_proc.set_ana_output_file("")

raw_viewer   = larlite.ClusterViewer()

my_proc.add_process(raw_viewer)

#raw_viewer.SetClusterProducer("fuzzycluster")
raw_viewer.SetClusterProducer("michel")
raw_viewer.SetHitProducer("cchit")

raw_viewer.SetZoomedView(True)

gStyle.SetOptStat(0)

#Check if handscan text file already exists:
does_outfile_exist = os.path.isfile(handscan_textfilename) 

#If outfile exists, read it in to a dictionary
scanned_events_dict = {}
if does_outfile_exist:
    print "Warning! Output file already exists. We will be appending to it. Delete if it you want to start fresh."
    for line in open(handscan_textfilename,'r'):
        columns = line.split(' ')
        myrun, mysubrun, myeventid, myidx = columns[0],columns[1],columns[2],columns[3]
        if (myrun, mysubrun, myeventid) not in scanned_events_dict.keys():
            scanned_events_dict[(myrun,mysubrun,myeventid)] = [myidx]
        else:
            scanned_events_dict[(myrun,mysubrun,myeventid)].append(myidx)
print scanned_events_dict
#Open the output text file
fout = open(handscan_textfilename,'a')

try: user_input_evt_no = int(input('Hit Enter to start on the first event, or type an int to jump to that event index.'))
except SyntaxError:
    user_input_evt_no = 0

while True:

    my_proc.process_event(user_input_evt_no)

    plane = 2
    for cindex in xrange(raw_viewer.ClusterCount(plane)):

        res = raw_viewer.DrawOneClusterGraphAndHits(plane,cindex)
        run, subrun, evtid, idx = int(res.runnumber), int(res.subrunnumber), int(res.eventid), int(res.index)
        if (run, subrun, evtid) in scanned_events_dict.keys():
              if idx in scanned_events_dict[(run, subrun, evtid)]:
                print "Whoops this cluster (run %d, subrun %d, eventid %d, index %d) has already been scanned. Skipping..." %(run, subrun, evtid, idx)
                continue

        while True:
            try: 
                good_or_bad = int(input('Type 0 if bad, 1 if good, 2 to skip cluster.'))
                if good_or_bad not in [0, 1, 2]:
                    print "You didn't input 0, 1, 2, or 3! Try again..."
                    continue
            except ValueError:
                print "You didn't input an integer."
                continue
            except NameError:
                print "You didn't input a number 1 (good) or zero (bad)."
                continue
            except SyntaxError:
                print "You didn't input anything at all. Try again."
                continue
            else:
                break
        if good_or_bad in [0, 1]:
            if (run, subrun, evtid) not in scanned_events_dict.keys():
                scanned_events_dict[(run, subrun, evtid)] = [idx]
            else:
                if idx not in scanned_events_dict[(run, subrun, evtid)]:
                    scanned_events_dict[(run, subrun, evtid)].append(idx)
            fout.write('%d %d %d %d %d\n'%(run, subrun, evtid, idx,good_or_bad))

    user_input_evt_no = user_input_evt_no + 1


