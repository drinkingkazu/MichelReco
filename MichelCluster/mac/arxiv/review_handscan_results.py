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
scanned_events_result_dict = {}
if does_outfile_exist:
    print "Warning! Output file already exists. We will be appending to it. Delete if it you want to start fresh."
    for line in open(handscan_textfilename,'r'):
        columns = line.split(' ')
        myrun, mysubrun, myeventid, myidx, myres = columns[0],columns[1],columns[2],columns[3],columns[4]
        myrun, mysubrun, myeventid, myidx, myres = int(myrun), int(mysubrun), int(myeventid), int(myidx), int(myres)
        if (myrun, mysubrun, myeventid, myidx) not in scanned_events_result_dict.keys():
            scanned_events_result_dict[(myrun,mysubrun,myeventid,myidx)] = myres
        else:
            print 'wtf you scanned this cluster twice? figure your shit out.'
            quit()


user_input_evt_no=0
while True:

    my_proc.process_event(user_input_evt_no)

    plane = 2
    for cindex in xrange(raw_viewer.ClusterCount(plane)):

        res = raw_viewer.DrawOneClusterGraphAndHits(plane,cindex)
        run, subrun, evtid, idx = int(res.runnumber), int(res.subrunnumber), int(res.eventid), int(res.index)
        if (run, subrun, evtid, idx) not in scanned_events_result_dict.keys():
            print "You never scanned this cluster index in this event. Skipping."
            continue
        else:      
            myres = scanned_events_result_dict[(run,subrun,evtid,idx)]
            good_or_bad = "GOOD" if myres == 1 else "BAD"
            print "This cluster was previously scanned to be %s" % good_or_bad
            if good_or_bad == "BAD":
                continue
            try: 
                meaningless = int(input('Hit Enter.'))
            except SyntaxError: 
                pass


    user_input_evt_no += 1



