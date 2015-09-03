import numpy as np
import ROOT
import matplotlib.pyplot as plt
from root_numpy import root2array, root2rec, tree2rec, array2root
import sys
import signal
import pandas as pd

if (len(sys.argv) < 2):
    print
    print "Incorrect Use"
    print "Please provide the path to the file containing the tree which stores info on the reco'ed Michel"
    print
    sys.exit(0)

# do you want to scan only throgh signal or background?
print
print "Do you want to scan only through Signal/Background events?"
print "[ \'s\' = signal. \'b\' = background. \'q\' = no thanks, scan everything]" 
print
# how many times have we tried getting the correct input from the user
attempts = 0
# parameter deciding if to scan signal (1) or background (0)
# value of 2 will scan everything
# set to -1 by default
whattoscan = -1
# text file where to search for signal/background handscan result
handscanresults = 'mac/michel_handscanning_results_v3_stepsonic.txt'
while (whattoscan == -1):
    ret = raw_input('Enter command...')
    if (ret == 's'):
        print
        print 'Great! Looking for signal events in %s'%handscanresults
        whattoscan = 1
    elif (ret == 'b'):
        print 'Great! Looking for background events in %s'%handscanresults
        whattoscan = 0
    elif (ret == 'q'):
        print 'Great! Looping through all events in %s'%handscanresults
        whattoscan = 2
    else:
        if (attempts >= 10):
            print 'Ok...we are done here...we will scan everything...'
            print
            whattoscan = 2
        else:
            print
            print 'You have failed at entering the requested input...try again...'
            print "Do you want to scan only through Signal/Background events?"
            print "[ \'s\' = signal. \'b\' = background. \'q\' = no thanks, scan everything]" 
            attempts += 1

scanDict = {}

if ((whattoscan == 0) or (whattoscan == 1)):
    fin = open(handscanresults)
    for line in fin:
        words = line.split()
        run    = int(words[0])
        subrun = int(words[1])
        event  = int(words[2])
        index  = int(words[3])
        what   = int(words[4])
        scanDict[(run,event,index)] = what

fin = sys.argv[-1]

#catch ctrl+c
def signal_handler(signal, frame):
    print "You pressed Ctrl+C. Exiting!"
    sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)

# interactive matplotlib
plt.ion()

plt.rcParams.update({'font.size': 14})

fig, axarr = plt.subplots(2,2,figsize=(13,13))

f = ROOT.TFile(fin)

t = f.Get('out_tree')



arr = tree2rec(t,branches=['_michel_clustered_charge','_michel_Z','_michel_X','_X','_Z','_q_v','_s_v','_chi_v','_t_dqds_v',
                           '_t_q_v','_mean_chi','_lowest_chi','_rms_chi','_boundary','_chi_at_boundary',
                           '_event','_run','_subrun','_clus_idx','_forward'])

all_hits = False
if (f.GetListOfKeys().Contains('_hit_tree')):
    all_hits = True
    #hit tree
    hit_tree = f.Get('_hit_tree')
    # get hit information
    hits = pd.DataFrame(tree2rec(hit_tree,branches=['_run','_subrun','_event','_q_v','_w_v','_t_v']))
    # fill a dictionary that goes from (run,subrun,event) -> data frame with hit info for that event
    hitDict = {}
    for n in xrange(len(hits)):
        run    = hits['_run'][n]
        subrun = hits['_subrun'][n]
        event  = hits['_event'][n]
        hitDict[(run,subrun,event)] = hits.query('_event == %i and _run == %i and _subrun == %i'%(event,run,subrun))


# define function that, given an array, returns a sub-array with values within a certain range
def getSubArray(xarr,xmin,xmax,yarr,ymin,ymax):
    if (len(xarr) != len(yarr)):
        return [],[]
    xret = []
    yret = []
    for n in xrange(len(xarr)):
        x = xarr[n]
        y = yarr[n]
        if ( (x > xmin) and (x < xmax) and (y > ymin) and (y < ymax) ):
            xret.append(x)
            yret.append(y)
    return xret,yret

# make list of bg entries
bgList = []
for x in xrange(len(arr)):
    evt = arr['_event'][x]
    run = arr['_run'][x]
    subrun = arr['_subrun'][x]
    idx = arr['_clus_idx'][x]
    if not (run,evt,idx) in scanDict:
        continue
    if (scanDict[(run,evt,idx)] < 2):
        bgList.append(x)

# number of event to scan
n = 0

for k in xrange(len(bgList)/4):

    print 'scanning TTree enrty: %i'%n

    # first make sure that if we are scanning through signal or background only
    # we are only looping through those events



    # clear the axes...
    axarr[0,0].cla()
    axarr[0,1].cla()
    axarr[1,0].cla()
    axarr[1,1].cla()
    fig.gca()

    #plt.tight_layout()

    evt = arr['_event'][bgList[n]]
    run = arr['_run'][bgList[n]]
    subrun = arr['_subrun'][bgList[n]]
    idx = arr['_clus_idx'][bgList[n]]

    clus_x = arr['_X'][bgList[n]]
    clus_z = arr['_Z'][bgList[n]]
    michel_x = arr['_michel_X'][bgList[n]]
    michel_z = arr['_michel_Z'][bgList[n]]
    boundary = arr['_boundary'][bgList[n]]
    forward = arr['_forward'][bgList[n]]
    dq    = arr['_q_v'][bgList[n]]
    ds    = arr['_s_v'][bgList[n]]
    tdq   = arr['_t_q_v'][bgList[n]]
    tdqds = arr['_t_dqds_v'][bgList[n]]

    # drawing michel cluster hit position
    #axarr[0,0].set_xlabel('x [cm]')
    #axarr[0,0].set_ylabel('z [cm]')
 
    # draw all hits in even
    # set bounds based on "clus" bounds
    xmax = np.amax(np.array(michel_x))
    xmin = np.amin(np.array(michel_x))
    zmax = np.amax(np.array(michel_z))
    zmin = np.amin(np.array(michel_z))
    dx = xmax-xmin
    dz = zmax-zmin
    xmax += 15
    xmin -= 15
    zmax += 15
    zmin -= 15
    if (all_hits):
        df = hitDict[(run,subrun,evt)]
        zpoints,xpoints = getSubArray(np.array(df['_w_v'])[0],zmin,zmax,np.array(df['_t_v'])[0],xmin,xmax)
        axarr[0,0].scatter(zpoints,xpoints,c='k',edgecolor='none',alpha=0.2,s=30)
    # draw full cluster
    axarr[0,0].scatter(clus_z,clus_x,c=dq,s=50,edgecolor='none')
    # draw michel cluster
    axarr[0,0].scatter(michel_z,michel_x,c='k',edgecolor='none',s=60)
    axarr[0,0].set_title('Evt: %i Run: %i Idx: %i'%(evt,run,idx))
    axarr[0,0].set_xlim([zmin,zmax])
    axarr[0,0].set_ylim([xmin,xmax])
    axarr[0,0].grid()


    ###################################################

    evt = arr['_event'][bgList[n+1]]
    run = arr['_run'][bgList[n+1]]
    subrun = arr['_subrun'][bgList[n+1]]
    idx = arr['_clus_idx'][bgList[n+1]]

    clus_x = arr['_X'][bgList[n+1]]
    clus_z = arr['_Z'][bgList[n+1]]
    michel_x = arr['_michel_X'][bgList[n+1]]
    michel_z = arr['_michel_Z'][bgList[n+1]]
    boundary = arr['_boundary'][bgList[n+1]]
    forward = arr['_forward'][bgList[n+1]]
    dq    = arr['_q_v'][bgList[n+1]]
    ds    = arr['_s_v'][bgList[n+1]]
    tdq   = arr['_t_q_v'][bgList[n+1]]
    tdqds = arr['_t_dqds_v'][bgList[n+1]]

    # drawing michel cluster hit position
    #axarr[0,1].set_xlabel('x [cm]')
    #axarr[0,1].set_ylabel('z [cm]')
 
    # draw all hits in even
    # set bounds based on "clus" bounds
    xmax = np.amax(np.array(michel_x))
    xmin = np.amin(np.array(michel_x))
    zmax = np.amax(np.array(michel_z))
    zmin = np.amin(np.array(michel_z))
    dx = xmax-xmin
    dz = zmax-zmin
    xmax += 15
    xmin -= 15
    zmax += 15
    zmin -= 15
    if (all_hits):
        df = hitDict[(run,subrun,evt)]
        zpoints,xpoints = getSubArray(np.array(df['_w_v'])[0],zmin,zmax,np.array(df['_t_v'])[0],xmin,xmax)
        axarr[0,1].scatter(zpoints,xpoints,c='k',edgecolor='none',alpha=0.2,s=30)
    # draw full cluster
    axarr[0,1].scatter(clus_z,clus_x,c=dq,s=50,edgecolor='none')
    # draw michel cluster
    axarr[0,1].scatter(michel_z,michel_x,c='k',edgecolor='none',s=60)
    axarr[0,1].set_title('Evt: %i Run: %i Idx: %i'%(evt,run,idx))
    axarr[0,1].set_xlim([zmin,zmax])
    axarr[0,1].set_ylim([xmin,xmax])
    axarr[0,1].grid()



    #################################

    evt = arr['_event'][bgList[n+2]]
    run = arr['_run'][bgList[n+2]]
    subrun = arr['_subrun'][bgList[n+2]]
    idx = arr['_clus_idx'][bgList[n+2]]

    clus_x = arr['_X'][bgList[n+2]]
    clus_z = arr['_Z'][bgList[n+2]]
    michel_x = arr['_michel_X'][bgList[n+2]]
    michel_z = arr['_michel_Z'][bgList[n+2]]
    boundary = arr['_boundary'][bgList[n+2]]
    forward = arr['_forward'][bgList[n+2]]
    dq    = arr['_q_v'][bgList[n+2]]
    ds    = arr['_s_v'][bgList[n+2]]
    tdq   = arr['_t_q_v'][bgList[n+2]]
    tdqds = arr['_t_dqds_v'][bgList[n+2]]

    # drawing michel cluster hit position
    #axarr[1,0].set_xlabel('x [cm]')
    #axarr[1,0].set_ylabel('z [cm]')
 
    # draw all hits in even
    # set bounds based on "clus" bounds
    xmax = np.amax(np.array(michel_x))
    xmin = np.amin(np.array(michel_x))
    zmax = np.amax(np.array(michel_z))
    zmin = np.amin(np.array(michel_z))
    dx = xmax-xmin
    dz = zmax-zmin
    xmax += 15
    xmin -= 15
    zmax += 15
    zmin -= 15
    if (all_hits):
        df = hitDict[(run,subrun,evt)]
        zpoints,xpoints = getSubArray(np.array(df['_w_v'])[0],zmin,zmax,np.array(df['_t_v'])[0],xmin,xmax)
        axarr[1,0].scatter(zpoints,xpoints,c='k',edgecolor='none',alpha=0.2,s=30)
    # draw full cluster
    axarr[1,0].scatter(clus_z,clus_x,c=dq,s=50,edgecolor='none')
    # draw michel cluster
    axarr[1,0].scatter(michel_z,michel_x,c='k',edgecolor='none',s=60)
    axarr[1,0].set_title('Evt: %i Run: %i Idx: %i'%(evt,run,idx))
    axarr[1,0].set_xlim([zmin,zmax])
    axarr[1,0].set_ylim([xmin,xmax])
    axarr[1,0].grid()




    ###################################

    evt = arr['_event'][bgList[n+3]]
    run = arr['_run'][bgList[n+3]]
    subrun = arr['_subrun'][bgList[n+3]]
    idx = arr['_clus_idx'][bgList[n+3]]

    clus_x = arr['_X'][bgList[n+3]]
    clus_z = arr['_Z'][bgList[n+3]]
    michel_x = arr['_michel_X'][bgList[n+3]]
    michel_z = arr['_michel_Z'][bgList[n+3]]
    boundary = arr['_boundary'][bgList[n+3]]
    forward = arr['_forward'][bgList[n+3]]
    dq    = arr['_q_v'][bgList[n+3]]
    ds    = arr['_s_v'][bgList[n+3]]
    tdq   = arr['_t_q_v'][bgList[n+3]]
    tdqds = arr['_t_dqds_v'][bgList[n+3]]

    # drawing michel cluster hit position
    #axarr[1,1].set_xlabel('x [cm]')
    #axarr[1,1].set_ylabel('z [cm]')
 
    # draw all hits in even
    # set bounds based on "clus" bounds
    xmax = np.amax(np.array(michel_x))
    xmin = np.amin(np.array(michel_x))
    zmax = np.amax(np.array(michel_z))
    zmin = np.amin(np.array(michel_z))
    dx = xmax-xmin
    dz = zmax-zmin
    xmax += 15
    xmin -= 15
    zmax += 15
    zmin -= 15
    if (all_hits):
        df = hitDict[(run,subrun,evt)]
        zpoints,xpoints = getSubArray(np.array(df['_w_v'])[0],zmin,zmax,np.array(df['_t_v'])[0],xmin,xmax)
        axarr[1,1].scatter(zpoints,xpoints,c='k',edgecolor='none',alpha=0.2,s=30)
    # draw full cluster
    axarr[1,1].scatter(clus_z,clus_x,c=dq,s=50,edgecolor='none')
    # draw michel cluster
    axarr[1,1].scatter(michel_z,michel_x,c='k',edgecolor='none',s=60)
    axarr[1,1].set_title('Evt: %i Run: %i Idx: %i'%(evt,run,idx))
    axarr[1,1].set_xlim([zmin,zmax])
    axarr[1,1].set_ylim([xmin,xmax])
    axarr[1,1].grid()



    fig.canvas
    
    plt.savefig('all_%03i.png'%n)

    fig.canvas.draw()

    usrinput = raw_input("Hit Enter: next evt  || int: go to event number ||  q: exit viewer\n")                        
    if ( usrinput == "q" ):                                                                                                                  
        sys.exit(0)
    elif ( usrinput.isdigit() == True ):
        n = int(usrinput)
    else:
        n += 4

sys.exit(0)
