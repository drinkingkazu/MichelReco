import sys, os

import ROOT
ROOT.gSystem.Load("libMichelReco_MichelClusterFmwk")
from larlite import larlite as fmwk
from ROOT import michel
from highStats_algo import PrepareMichelAlgo

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file

#PATH = '/uboone/data/users/davidc1/Michel/MCcosmic/translate/data/4527415_'
PATH = '/uboone/data/users/davidc1/Michel/Paddles_Run3336/'

for x in xrange(0,50):

    #fname_clusters = 'larlite_rawclusters_%04i.root'%x
    fname_clusters = 'larlite_clusters_%04i.root'%x
    
    file_path_clusters = '%s%s'%(PATH,fname_clusters)

    print 'file is %s'%file_path_clusters

    if (os.path.isfile(file_path_clusters) == False):
        continue

    # use rawcluster files
    #if (file_path.find('raw') < 0):
    #    continue
    # use reco'd cluster files
    if (file_path_clusters.find('larlite_cluster') < 0):
        continue

    #### ADDING GENERIC DATA ####
    my_proc.add_input_file(file_path_clusters)

    fname_reco = 'larlite_opdata_%04i.root'%x
    
    file_path_reco = '%s%s'%(PATH,fname_reco)

    print 'file is %s'%file_path_reco

    if (os.path.isfile(file_path_reco) == False):
        continue

    if (file_path_reco.find('larlite_opdata') < 0):
        continue

    my_proc.add_input_file(file_path_reco)


    
    my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
    #my_proc.set_io_mode(fmwk.storage_manager.kREAD)

    my_proc.set_ana_output_file("%s/michel_tree_%04i.root"%(PATH,x))
    outfile = "%s/michel_clusters_%04i.root"%(PATH,x)
    my_proc.set_output_file( outfile )
    
    print 'output file is %s'%outfile

    #########################
    # Michel reco driver code
    my_unit = fmwk.MichelRecoDriver()
    my_unit.SetClusterProducer("fuzzycluster")
    #my_unit.SetClusterProducer("rawcluster")
    #my_unit.SetClusterProducer("linecluster")
    #my_unit.saveMichelClusters(True)
    my_unit.saveMichelClusters(True)
    #my_unit.SetClusterProducer("linecluster")
    my_unit.SetUseMC(False)
    my_unit.SetEField(0.27)
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

    ##############################
    # Attach algorithm for merging
    #mgr.AddMergingAlgo(michel.EdgeMerger())


    algoList = PrepareMichelAlgo()

    for algo in algoList:
        mgr.AddAlgo(algo)

    # Attach ana unit
    mgr.AddAna(michel.CosmicAna())

    # add process
    my_proc.add_process(my_unit)

    my_proc.set_data_to_write(fmwk.data.kCluster,'michel')
    my_proc.set_data_to_write(fmwk.data.kCluster,'muon')
    my_proc.set_data_to_write(fmwk.data.kAssociation,'michel')
    my_proc.set_data_to_write(fmwk.data.kAssociation,'muon')

    my_proc.enable_event_alignment(False)
    print
    print  "Finished configuring ana_processor. Start event loop!"
    print

    my_proc.run()

    # done!
    print
    print "Finished running ana_processor event loop!"
    print

sys.exit(0)
