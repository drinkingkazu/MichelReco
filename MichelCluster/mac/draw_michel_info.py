import numpy as np
import ROOT
import matplotlib
#matplotlib.use('TkAgg')
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
#handscanresults = 'mac/michel_handscanning_results_v3_stepsonic.txt'
handscanresults = 'mac/michel_handscanning_results.txt'
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
        print what

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
scale2 = axarr[1,1].twinx()

f = ROOT.TFile(fin)

t   = f.Get( 'out_tree' )


arr = tree2rec(t,branches=['_michel_clustered_charge',
                           '_michel_Z','_michel_X','_X',
                           '_Z','_q_v','_s_v','_chi_v','_t_dqds_v',
                           '_t_q_v','_mean_chi','_lowest_chi',
                           '_rms_chi','_boundary','_chi_at_boundary',
                           '_event','_run','_subrun','_clus_idx','_forward'])

ismc = False;
if f.GetListOfKeys().Contains('_mc_tree'):
    ismc = True

all_hits = False
if (f.GetListOfKeys().Contains('_hit_tree')):
    all_hits = True
    #hit tree
    hit_tree    = f.Get('_hit_tree')


    # get hit information
    hits       = pd.DataFrame(tree2rec(hit_tree   ,branches=['_run','_subrun','_event','_q_v','_w_v','_t_v','_p_v']))

    hitDict = {}
        
    hit_mc_tree = ROOT.TTree()
    mchit_info  = pd.DataFrame()

    if ismc :
        hit_mc_tree = f.Get( '_mc_tree'  )
        mchit_info  = pd.DataFrame(tree2rec(hit_mc_tree,branches=['_run','_subrun','_event','_michel_hit_frac']))
        hits        = pd.concat([hits,mchit_info._michel_hit_frac],axis=1)

    # fill a dictionary that goes from (run,subrun,event) -> data frame with hit info for that event
        
    for n in xrange(len(hits)):
        run    = hits['_run'][n]
        subrun = hits['_subrun'][n]
        event  = hits['_event'][n]
        hitDict[(run,subrun,event)] = hits.query('_event == %i and _run == %i and _subrun == %i'%(event,run,subrun))


# define function that, given an array, returns a sub-array with values within a certain range
def getSubArray(xarr,xmin,xmax,yarr,ymin,ymax,parr):
    if (len(xarr) != len(yarr)):
        return [],[]
    xret = []
    yret = []
    idx  = []
    for n in xrange(len(xarr)):
        x = xarr[n]
        y = yarr[n]
        p = parr[n]
        if ( (x > xmin) and (x < xmax) and (y > ymin) and (y < ymax) and p == 2):
            xret.append(x)
            yret.append(y)
            idx.append(n)
            
    #return xret, yret, idx
    return idx

# number of event to scan
n = 0
print len(arr)
while ( n < len(arr) ):

    print 'scanning TTree enrty: %i'%n

    # first make sure that if we are scanning through signal or background only
    # we are only looping through those events

    evt    = arr['_event'][n]
    run    = arr['_run'][n]
    subrun = arr['_subrun'][n]
    idx    = arr['_clus_idx'][n]



    if (whattoscan == 0):
        # if it cannot even find entry, continue
        if not (run,evt,idx) in scanDict:
            n += 1
            continue
        if (scanDict[(run,evt,idx)] != 0):
            n += 1
            continue
    if (whattoscan == 1):
        # if it cannot even find entry, continue
        if not (run,evt,idx) in scanDict:
            n += 1
            continue
        if (scanDict[(run,evt,idx)] != 1):
            n += 1
            continue

    # clear the axes...
    axarr[0,0].cla()
    axarr[0,1].cla()
    axarr[1,0].cla()
    axarr[1,1].cla()
    scale2.cla()
    fig.gca()

    clus_x = arr['_X'][n]
    clus_z = arr['_Z'][n]
    michel_x = arr['_michel_X'][n]
    michel_z = arr['_michel_Z'][n]

    dq    = arr['_q_v'][n]
    ds    = arr['_s_v'][n]
    tdq   = arr['_t_q_v'][n]
    tdqds = arr['_t_dqds_v'][n]
    
    boundary = arr['_boundary'][n]
    forward = arr['_forward'][n]

    Qtot = arr['_michel_clustered_charge'][n]

    # drawing chi vector vs. distance from beginning of cluster
    chi_v = arr['_chi_v'][n]
    chi_rms = arr['_rms_chi'][n]
    chi_mean = arr['_mean_chi'][n]
    chi_low = arr['_lowest_chi'][n]
    chi_boundary = arr['_chi_at_boundary'][n]
    # calculate a chi for the michel and one for the muon
    chi_forward = 0
    chi_backward   = 0
    for x in xrange(boundary,len(chi_v)):
        chi_forward += chi_v[x]
    if (chi_forward != 0):
        chi_forward /= float((len(chi_v)-boundary))
    for x in xrange(boundary):
        chi_backward += chi_v[x]
    if (chi_backward != 0):
        chi_backward /= float(boundary)

    # vector to hold charge for michel
    if (forward):
        chi_muon = chi_backward;
        chi_michel = chi_forward;
    else:
        chi_muon = chi_forward;
        chi_michel = chi_backward;

    axarr[1,0].set_xlabel('S [cm]')
    axarr[1,0].set_ylabel('Chi per Hit')
    #l = 'rms: %.02f mean: %.02f\nlow: %.02f bound: %.02f'%(chi_rms,chi_mean,chi_low,chi_boundary))
    l = 'chi michel: %.02f \nchi muon: %.02f'%(chi_michel,chi_muon)
    axarr[1,0].plot(ds,chi_v,'bo--',label=l)
    axarr[1,0].legend(loc=4)
    axarr[1,0].axvline(ds[boundary],lw=3,color='r')
    axarr[1,0].set_ylim([0,1])
    if (forward):
        axarr[1,0].axvspan(ds[boundary],ds[-1],color='r',alpha=0.2)
    else:
        axarr[1,0].axvspan(ds[0],ds[boundary],color='r',alpha=0.2)
    axarr[1,0].grid()

    # drawing michel cluster hit position
    axarr[0,0].set_xlabel('x [cm]')
    axarr[0,0].set_ylabel('z [cm]')
 
    # draw all hits in even
    # set bounds based on "clus" bounds
    if ismc:
        xmax = np.amax(np.array(clus_x))
        xmin = np.amin(np.array(clus_x))
        zmax = np.amax(np.array(clus_z))
        zmin = np.amin(np.array(clus_z))
    else:
        xmax = np.amax(np.array(michel_x))
        xmin = np.amin(np.array(michel_x))
        zmax = np.amax(np.array(michel_z))
        zmin = np.amin(np.array(michel_z))
    
    dx = xmax-xmin
    dz = zmax-zmin

    if ismc:
        xmax += dx/2
        xmin -= dx/2
        zmax += dz/2
        zmin -= dz/2
    else:
        xmax += 15
        xmin -= 15
        zmax += 15
        zmin -= 15
    print evt
    print subrun
    print run
    
    if (all_hits):
        df   = hitDict[(run,subrun,evt)]
        idxx = getSubArray(np.array(df['_w_v'])[0],
                          zmin,zmax,
                          np.array(df['_t_v'])[0],
                          xmin,xmax,
                          df['_p_v'].values[0])
        
        zpoints = df['_w_v'].values[0][idxx]
        xpoints = df['_t_v'].values[0][idxx]

        weights = 30

        cc = 'k'
        
        if ismc:
            weights = 200*df['_michel_hit_frac'].values[0][idxx]
            cc='purple'
        
        axarr[0,0].scatter(zpoints,
                           xpoints,
                           c=cc,
                           edgecolor='none',
                           alpha=0.75,
                           s=weights)
        
        

        
    # draw full cluster
    axarr[0,0].scatter(clus_z,clus_x,c=dq,s=50,edgecolor='none')
    # draw michel cluster
    axarr[0,0].scatter(michel_z,michel_x,c='k',edgecolor='none',s=60)

    axarr[0,0].set_title("Evt: %i Run: %i Idx: %i" % (evt,run,idx) )
    axarr[0,0].set_xlim([zmin,zmax])
    axarr[0,0].set_ylim([xmin,xmax])
    axarr[0,0].grid()


    # drawing charge vs. distance from beginning of cluster
    axarr[0,1].set_xlabel('S [cm]')
    axarr[0,1].set_ylabel('Q')
    axarr[0,1].plot(ds,dq,'bo--')
    axarr[0,1].axvline(ds[boundary],lw=3,color='r')
    if (forward):
        axarr[0,1].axvspan(ds[boundary],ds[-1],color='r',alpha=0.2)
    else:
        axarr[0,1].axvspan(ds[0],ds[boundary],color='r',alpha=0.2)
    axarr[0,1].grid()



    # drawing truncated Q vs. distance from beginning of cluster
    # also truncated dQds vs. distance from beginning of cluster
    axarr[1,1].set_xlabel('S [cm]')
    axarr[1,1].set_ylabel('Truncated Q',color='b')
    axarr[1,1].plot(ds,tdq,'bo--')
    axarr[1,1].set_ylim(axarr[1,1].get_ylim())
    for ti in axarr[1,1].get_yticklabels():
        ti.set_color('b')

    scale2.plot(ds,tdqds,'mo--')
    scale2.set_ylabel('Truncated dQ/ds',color='m')
    for ti in scale2.get_yticklabels():
        ti.set_color('m')
    axarr[1,1].axvline(ds[boundary],lw=3,color='r')
    if (forward):
        axarr[1,1].axvspan(ds[boundary],ds[-1],color='r',alpha=0.2)
    else:
        axarr[1,1].axvspan(ds[0],ds[boundary],color='r',alpha=0.2)
    axarr[1,1].grid()

    #plt.show()
    fig.canvas
    fig.canvas.draw()

    usrinput = raw_input("Hit Enter: next evt  || int: go to event number ||  q: exit viewer\n")                        
    if ( usrinput == "q" ):                                                                                                                  
        sys.exit(0)
    elif ( usrinput.isdigit() == True ):
        n = int(usrinput)
    else:
        n += 1

sys.exit(0)
