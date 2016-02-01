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


f = ROOT.TFile(sys.argv[-1])
t   = f.Get( 'out_tree' )


arr = tree2rec(t,branches=['_michel_clustered_charge',
                           '_michel_n_hits',
                           'michel_start_X','michel_start_Z',
                           '_michel_Z','_michel_X','_X',
                           '_Z','_q_v','_s_v','_chi_v','_t_dqds_v',
                           '_t_q_v','_mean_chi','_lowest_chi',
                           '_rms_chi','_boundary','_chi_at_boundary',
                           '_event','_run','_subrun','_clus_idx','_forward'])


icarus_f = open('/home/david/Desktop/icarus_data.csv','r')

icarus_energies = []
icarus_numbers  = []
icarus_errors   = []
icarus_area = 0

e_prev = 0

for line in icarus_f:
    words = line.split()
    e = float(words[0])
    n = float(words[1])
    err= np.sqrt(float(n))
    icarus_energies.append(e)
    icarus_numbers.append(n)
    icarus_errors.append(err)
    icarus_area += n*(e-e_prev)
    e_prev = e

# normalize
icarus_energies = np.array(icarus_energies)
icarus_numbers = np.array(icarus_numbers)
icarus_errors = np.array(icarus_errors)
icarus_numbers /= float(icarus_area)
icarus_errors /= float(icarus_area)

michel_energies = []

# number of event to scan
n = 0

while ( n < len(arr) ):

    # first make sure that if we are scanning through signal or background only
    # we are only looping through those events
    evt    = arr['_event'][n]
    run    = arr['_run'][n]
    subrun = arr['_subrun'][n]
    idx    = arr['_clus_idx'][n]

    dq    = np.array(arr['_q_v'][n])

    maxdq = np.amax(dq)
    
    Qtot  = np.array(arr['_michel_clustered_charge'][n])
    Nhits = np.array(arr['_michel_n_hits'][n])

    n += 1

    #print 'Michel Charge [ADCs] : %.02f'%Qtot
    #print 'Max Q hit [ADCs]     : %.02f'%maxdq
    #print 'Michel Energy [MeV]  : %.02f'%EMeV
    #print 'Michel Tot Hits      : %i'%Nhits
    #print 'Charge / Hit         : %.02f'%(Qtot/Nhits)

    Qmichel = Qtot-maxdq
    Emichel = Qmichel * 0.00263 * 2 * 1.15

    if (Qtot/Nhits > 400):
        continue



    michel_energies.append(Emichel)



fig = plt.figure(figsize=(10,10))

emin  = 0
nbins = 40
emax  = 100
ebin = (emax-emin)/float(nbins)
emin -= ebin/2.
emax += ebin/2.
Ebins = np.linspace(emin,emax,nbins)
entries, bins = np.histogram(michel_energies,Ebins)
entries = entries.astype(float)
bin_centers = 0.5*(bins[1:]+bins[:-1])

# calculate normalization factor from integral
tot_michels = np.sum(entries)
bin_width = bin_centers[1]-bin_centers[0]
area = tot_michels * bin_width
print 'found %i michels'%tot_michels

y_err = np.sqrt(entries)
entries /= float(area)
y_err /= float(area)

plt.errorbar(bin_centers,entries,yerr=y_err,fmt='o',color='b',label='MicroBooNE',markersize=8)
#entries,bins, etc = plt.hist(michel_energies,Ebins,histtype='stepfilled',color='b',normed=True)
plt.grid()
plt.xlabel('Energy [MeV]')
plt.ylabel('Relative Count [Area Normalized]')
plt.title('Reconstructed Michel Energy Spectrum')
plt.xlim([0,80])
plt.errorbar(icarus_energies,icarus_numbers,yerr=icarus_errors,fmt='o',color='r',label='ICARUS T600',markersize=8)
plt.legend(loc=1,numpoints=1)
plt.show()
