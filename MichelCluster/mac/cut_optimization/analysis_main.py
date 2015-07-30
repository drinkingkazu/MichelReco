from analysis_helper import *
from analysis_loader import *


print 'Reading in signal...'
SIGNAL     = get_data_frame("/Users/vgenty/Dropbox/ana3_michel.root","out_tree")

print 'Reading in background...'
BACKGROUND = get_data_frame("/Users/vgenty/Dropbox/ana4_muon.root","out_tree")

############Apply any pre-cuts############
print 'Applying precuts...'
print ''
print 'Signal'

S = SIGNAL.query('_largest_cluster_charge > 0')
S = S.query('_number_of_clusters >= 1')
S = S.query('_n_hits_in_largest_cluster_michel >= 10')
# S = S.query('lowest_chi < 0.25')
# S = S.query('chi_at_boundary < 0.4')
# S = S.query('mean_chi > 0.8')

print 'Background'
B = BACKGROUND.query('_largest_cluster_charge > 0')
B = B.query('_number_of_clusters >= 1')
B = B.query('_n_hits_in_largest_cluster_michel >= 10')
# B = B.query('lowest_chi < 0.25')
# B = B.query('chi_at_boundary < 0.4')
# B = B.query('mean_chi > 0.8')

print ''
##########################################

print 'Creating data...'
data = DataHolder(4828,  #nmuons
                  61692, #nmichels
                  150,   #bins
                  0,     #xlow
                  150,   #xhigh
                  S, 
                  B,
                  '_largest_cluster_charge')


print "Creating cuts..."
#.....! 2 Above
# cut1Name = '_n_hits_in_largest_cluster_michel'
# cut2Name = '_largest_cluster_charge'
# cut1Vars  = np.arange(0,50, 1)
# cut2Vars  = np.arange(0,100,2)


#Good! 1st below, 1st above..
# cut1Name = 'lowest_chi'
# cut2Name = '_n_hits_in_largest_cluster'
# cut1Vars  = np.arange(0,1,0.05)
# cut2Vars  = np.arange(0,100,5)

#Good 1st below, 1st above..
# cut1Name = 'lowest_chi'
# cut2Name = 'mean_chi'
# cut1Vars  = np.arange(0,1,0.05)
# cut2Vars  = np.arange(0,1,0.05)

#Ok... 1st below, 1st above..
# cut1Name = 'chi_at_boundary'
# cut2Name = 'mean_chi'
# cut1Vars  = np.arange(0,1,0.05)
# cut2Vars  = np.arange(0,1,0.05)

# #....... 2 below
cut1Name = 'chi_at_boundary'
cut2Name = 'lowest_chi'
cut1Vars  = np.arange(0,1,0.01)
cut2Vars  = np.arange(0,1,0.01)

# cut1Name = 'lowest_chi'
# cut2Name = '_n_hits_in_largest_cluster_michel'
# cut1Vars  = np.arange(0,1,0.01)
# cut2Vars  = np.arange(0,50, .5)



#.....! 2 Above
# cut1Name = 'mean_chi'
# cut2Name = '_n_hits_in_largest_cluster_michel'
# cut1Vars = np.arange(0,1,0.05)
# cut2Vars = np.arange(0,50, 1)



print "Computing sensitivity between... %s and %s" % (cut1Name, cut2Name)
# Ldata  = [[calc_sig(data.get_histogram_2cuts_S_2_above(cut1Name,
#                                                        cut2Name,
#                                                        cut1,
#                                                        cut2,
#                                                        150,
#                                                        150),
#                     data.get_histogram_2cuts_B_2_above(cut1Name,
#                                                        cut2Name,
#                                                        cut1,
#                                                        cut2,
#                                                        150,
#                                                        150))
#            for cut1 in cut1Vars] for cut2 in cut2Vars]
# Ldata  = [[calc_sig(data.get_histogram_2cuts_S_1_above_1_below(cut1Name,
#                                                                cut2Name,
#                                                                cut1,
#                                                                cut2,
#                                                                150,
#                                                                150),
#                     data.get_histogram_2cuts_B_1_above_1_below(cut1Name,
#                                                                cut2Name,
#                                                                cut1,
#                                                                cut2,
#                                                                150,
#                                                                150))
#            for cut1 in cut1Vars] for cut2 in cut2Vars]
Ldata  = [[calc_sig(data.get_histogram_2cuts_S_2_below(cut1Name,
                                                       cut2Name,
                                                       cut1,
                                                       cut2,
                                                       150,
                                                       150),
                    data.get_histogram_2cuts_B_2_below(cut1Name,
                                                       cut2Name,
                                                       cut1,
                                                       cut2,
                                                       150,
                                                       150))
           for cut1 in cut1Vars] for cut2 in cut2Vars]


plt.figure(figsize=(12,6))

plt.contourf(cut1Vars,cut2Vars,Ldata - data.total_sig,200) #200 is bin width ?
plt.xlabel(cut1Name,fontsize=20)
plt.ylabel(cut2Name,fontsize=20)

plt.colorbar()
plt.tick_params(labelsize=20)
plt.show()
