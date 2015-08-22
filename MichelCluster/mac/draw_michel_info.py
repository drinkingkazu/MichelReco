import numpy as np
import ROOT
import matplotlib.pyplot as plt
from root_numpy import root2array, root2rec, tree2rec, array2root
import sys

if (len(sys.argv) < 2):
    print
    print "Incorrect Use"
    print "Please provide the path to the file containing the tree which stores info on the reco'ed Michel"
    print
    sys.exit(0)


fin = sys.argv[-1]

# interactive matplotlib
plt.ion()

plt.rcParams.update({'font.size': 14})

fig, axarr = plt.subplots(2,2,figsize=(13,13))
scale2 = axarr[1,1].twinx()

f = ROOT.TFile(fin)

t = f.Get('out_tree')

arr = tree2rec(t,branches=['_michel_clustered_charge','_michel_Z','_michel_X','_X','_Z','_q_v','_s_v','_chi_v','_t_dqds_v',
                           '_t_q_v','_mean_chi','_lowest_chi','_rms_chi','_boundary','_chi_at_boundary',
                           '_event','_run','_subrun','_clus_idx'])

for n in xrange(len(arr)):

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

    dq = arr['_q_v'][n]
    ds = arr['_s_v'][n]
    tdq = arr['_t_q_v'][n]
    tdqds = arr['_t_dqds_v'][n]
    
    boundary = arr['_boundary'][n]

    Qtot = arr['_michel_clustered_charge'][n]

    evt = arr['_event'][n]
    run = arr['_run'][n]
    idx = arr['_clus_idx'][n]

    # drawing michel cluster hit position
    axarr[0,0].set_xlabel('x [cm]')
    axarr[0,0].set_ylabel('z [cm]')
    axarr[0,0].scatter(clus_x,clus_z,c='k',s=10)
    axarr[0,0].scatter(michel_x,michel_z,c='r',edgecolor='none',s=30)
    axarr[0,0].set_title('Evt: %i Run: %i Idx: %i'%(evt,run,idx))
    axarr[0,0].grid()


    # drawing charge vs. distance from beginning of cluster
    axarr[0,1].set_xlabel('S [cm]')
    axarr[0,1].set_ylabel('Q')
    axarr[0,1].plot(ds,dq,'bo--')
    axarr[0,1].axvline(ds[boundary],lw=3,color='r')
    axarr[0,1].grid()

    # drawing chi vector vs. distance from beginning of cluster
    chi_v = arr['_chi_v'][n]
    chi_rms = arr['_rms_chi'][n]
    chi_mean = arr['_mean_chi'][n]
    chi_low = arr['_lowest_chi'][n]
    chi_boundary = arr['_chi_at_boundary'][n]
    axarr[1,0].set_xlabel('S [cm]')
    axarr[1,0].set_ylabel('Chi per Hit')
    axarr[1,0].plot(ds,chi_v,'bo--',label='rms: %.02f mean: %.02f\nlow: %.02f bound: %.02f'%(chi_rms,chi_mean,chi_low,chi_boundary))
    axarr[1,0].legend(loc=4)
    axarr[1,0].axvline(ds[boundary],lw=3,color='r')
    axarr[1,0].grid()

    # drawing truncated Q vs. distance from beginning of cluster
    # also truncated dQds vs. distance from beginning of cluster
    axarr[1,1].set_xlabel('S [cm]')
    axarr[1,1].set_ylabel('Truncated Q',color='b')
    axarr[1,1].plot(ds,tdq,'bo--')
    for ti in axarr[1,1].get_yticklabels():
        ti.set_color('b')

    scale2.plot(ds,tdqds,'mo--')
    scale2.set_ylabel('Truncated dQ/ds',color='m')
    for ti in scale2.get_yticklabels():
        ti.set_color('m')
    axarr[1,1].axvline(ds[boundary],lw=3,color='r')
    axarr[1,1].grid()

    plt.show()

    usrinput = raw_input("Hit Enter: next evt  ||  q: exit viewer\n")                                                                        
    if ( usrinput == "q" ):                                                                                                                  
        sys.exit(0)

sys.exit(0)
