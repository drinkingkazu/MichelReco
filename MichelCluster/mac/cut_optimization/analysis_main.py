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
S = S.query('_n_hits_in_largest_cluster_michel >= 1')
# S = S.query('lowest_chi < 0.25')
# S = S.query('chi_at_boundary < 0.4')
# S = S.query('mean_chi > 0.8')

print 'Background'
B = BACKGROUND.query('_largest_cluster_charge > 0')
B = B.query('_number_of_clusters >= 1')
B = B.query('_n_hits_in_largest_cluster_michel >= 1')
# B = B.query('lowest_chi < 0.25')
# B = B.query('chi_at_boundary < 0.4')
# B = B.query('mean_chi > 0.8')

print ''
##########################################

print 'Creating data...'
data = DataHolder(4828,  #nmichels
                  61692, #nmuons
                  150,   #bins
                  0,     #xlow
                  150,   #xhigh
                  S, 
                  B,
                  '_largest_cluster_charge')




print "Creating cuts..."
cuts = {'_n_hits_in_largest_cluster_michel' : np.arange(0,50,1),
        '_largest_cluster_charge'           : np.arange(0,100,2),
        'lowest_chi'                        : np.arange(0,1,0.05),
        'mean_chi'                          : np.arange(0,1,0.05),
        'chi_at_boundary'                   : np.arange(0,1,0.05)
        }

print "Choose from possible cuts..."
for key in cuts:
    print key


cut1Name = ''
cut2Name = ''
cut1Vars = ''
cut2Vars = ''

while True:
    try:
        cut1Name = str(raw_input('Cut1 Name: '))
        cut2Name = str(raw_input('Cut2 Name: '))
        cut1Vars = cuts[cut1Name]
        cut2Vars = cuts[cut2Name]
    except KeyError:
        print 'Chosen key did not exist in list of'
        print 'available cuts. Try again.'
        continue

    break
    

# cut1Name = 'chi_at_boundary'
# cut2Name = 'lowest_chi'
# cut3Name = 'mean_chi'
# cut1Vars  = np.arange(0,1,0.02)
# cut2Vars  = np.arange(0,1,0.02)
# cut3Vars  = np.arange(0,1,0.02)




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
#                                                                  cut2Name,
#                                                                  cut1,
#                                                                  cut2,
#                                                                  150,
#                                                                  150),
#                       data.get_histogram_2cuts_B_1_above_1_below(cut1Name,
#                                                                  cut2Name,
#                                                                  cut1,
#                                                                  cut2,
#                                                                  150,
#                                                                  150))
#            for cut1 in cut1Vars] for cut2 in cut2Vars]
# Ldata  = [[calc_sig(data.get_histogram_2cuts_S_2_below(cut1Name,
#                                                        cut2Name,
#                                                        cut1,
#                                                        cut2,
#                                                        150,
#                                                        150),
#                     data.get_histogram_2cuts_B_2_below(cut1Name,
#                                                        cut2Name,
#                                                        cut1,
#                                                        cut2,
#                                                        150,
#                                                        150))
#            for cut1 in cut1Vars] for cut2 in cut2Vars]

#    def get_histogram_3cuts_S_2_below_1_above(self,
Ldata  = [[[calc_sig(data.get_histogram_3cuts_S_2_below_1_above(cut1Name,
                                                                cut2Name,
                                                                cut3Name,
                                                                cut1,
                                                                cut2,
                                                                cut3,
                                                                150,
                                                                150),
                     data.get_histogram_3cuts_B_2_below_1_above(cut1Name,
                                                                cut2Name,
                                                                cut3Name,
                                                                cut1,
                                                                cut2,
                                                                cut3,
                                                                150,
                                                                150))
            for cut1 in cut1Vars] for cut2 in cut2Vars] for cut3 in cut3Vars]



Ldata = np.array(Ldata)
Ldata -= data.total_sig;

from numpy import unravel_index
lowest = unravel_index(Ldata.argmax(),Ldata.shape)

# print 'The maxindicies in this histo are: (%d,%d)' % (lowest[1],
#                                                       lowest[0])

# print 'cut1Vars max: ',cut1Vars[lowest[1]]
# print 'cut2Vars max: ',cut2Vars[lowest[0]]

# plt.figure(figsize=(12,6))

# plt.contourf(cut1Vars,cut2Vars,Ldata,200) #200 is bin width ?
# plt.xlabel(cut1Name,fontsize=20)
# plt.ylabel(cut2Name,fontsize=20)

# plt.colorbar()
# plt.tick_params(labelsize=20)
# plt.show()


