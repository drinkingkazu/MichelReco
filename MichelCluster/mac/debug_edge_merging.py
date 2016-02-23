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
from highStats_algo import PrepareMichelAlgo

plt.rcParams.update({'font.size':16})

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
my_unit.SetClusterProducer("pandoraCosmic")
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

merge_algo = michel.EdgeMerger()
merge_algo.SetEdgeDistance(2)
mgr.AddMergingAlgo(merge_algo)

#algoList = PrepareMichelAlgo()
#for algo in algoList:
#    mgr.AddAlgo(algo)

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

plt.ion()
fig, ax = plt.subplots(figsize=(10,10))

while my_proc.process_event():

    # get input clusters
    input_clusters = mgr.GetInputClusters()

    # for every cluster pair, run the EdgeMerger algorithm
    # and visualize each merge
    for i in xrange(input_clusters.size()):
        for j in xrange(i+1,input_clusters.size()):

            ax.cla()
            fig.gca()
            
            cl1 = input_clusters.at(i)
            cl2 = input_clusters.at(j)

            merge_result = merge_algo.Merge(cl1,cl2)

            if (merge_result == False):
                continue
            
            hits1 = cl1._all_hits;
            hits2 = cl2._all_hits;

            print 'cluster 1 : linear / all hits = %.02f'%(float(len(cl1._hits))/len(cl1._all_hits))
            print 'cluster 2 : linear / all hits = %.02f'%(float(len(cl2._hits))/len(cl2._all_hits))

            times_1 = []
            wires_1 = []
            for h in hits1:
                times_1.append(h._t)
                wires_1.append(h._w)

            times_2 = []
            wires_2 = []
            for h in hits2:
                times_2.append(h._t)
                wires_2.append(h._w)

            plt.plot(wires_1,times_1,'ro',markeredgewidth=None,markersize=10)
            plt.plot(cl1._start._w,cl1._start._t,'*',color='r',markersize=15)
            plt.plot(cl1._end._w,cl1._end._t,'D',color='r',markersize=15)
            plt.plot(wires_2,times_2,'bo',markeredgewidth=None,markersize=10)
            plt.plot(cl2._start._w,cl2._start._t,'*',color='b',markersize=15)
            plt.plot(cl2._end._w,cl2._end._t,'D',color='b',markersize=15)
            plt.xlabel('Wire Coordinate [cm]')
            plt.ylabel('Time Coordinate [cm]')

            plt.grid()
            if (merge_result == True):
                plt.title('Merged Two Input Clusters')
                print 'merged clusters of size %i, %i'%(hits1.size(),hits2.size())
            else:
                plt.title('Not Merged!')
            #plt.show()

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
