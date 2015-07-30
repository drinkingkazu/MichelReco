#vgenty

import matplotlib.pyplot as plt
import numpy as np

from analysis_constants import *



def calc_sig(S,B):
    sig = 0
    for i in xrange(len(S)): 
        if (S[i] + B[i]) == 0:
           continue
        sig += S[i]/np.sqrt(S[i]+B[i])
    return sig




class DataHolder:
    def __init__(self,n_michels,n_muons,bins,xlow,xhigh,michelDF,muonDF,parameter):
        self.nmichels = n_michels;
        self.nmuons   = n_muons;
        self.bins     = bins;
        self.xhigh    = xhigh;
        self.xlow     = xlow;
        
        #S = signal
        #B = background
        
        self.muonDF    = muonDF;
        self.michelDF  = michelDF;
        self.parameter = parameter
        
        (self.nS, self.binsS, patchesS) = plt.hist(np.array(michelDF[parameter])*QtoE,
                                                   bins=self.bins,
                                                   range=(self.xlow,self.xhigh))
        
        (self.nB, self.binsB, patchesB) = plt.hist(np.array(muonDF[parameter])*QtoE,
                                                   bins=self.bins,
                                                   range=(self.xlow,self.xhigh))
        
        print 'Calculating total sensitivity'
        self.total_sig = calc_sig(self.nS,[0 for i in xrange(len(self.nS))])
        print 'Total sentivity is ',self.total_sig
        
    def get_histogram_2cuts_B_2_above(self,
                                      cut1Name,
                                      cut2Name,
                                      cut1,
                                      cut2,
                                      xhigh,
                                      binz) :
        
        return np.histogram(self.muonDF[(self.muonDF[cut1Name] >= cut1 ) & (self.muonDF[cut2Name]*QtoE >= cut2)][self.parameter]*QtoE,
                            bins=binz,
                            range=(0,xhigh))[0]

    def get_histogram_2cuts_S_2_above(self,
                                      cut1Name,
                                      cut2Name,
                                      cut1,
                                      cut2,
                                      xhigh,
                                      binz) :
        
        return np.histogram(self.michelDF[(self.michelDF[cut1Name] >= cut1 ) & (self.michelDF[cut2Name]*QtoE >= cut2)][self.parameter]*QtoE,
                            bins=binz,
                            range=(0,xhigh))[0]


    def get_histogram_2cuts_B_1_above_1_below(self,
                                              cut1Name,
                                              cut2Name,
                                              cut1,
                                              cut2,
                                              xhigh,
                                              binz) :
        
        return np.histogram(self.muonDF[(self.muonDF[cut1Name] <= cut1 ) & (self.muonDF[cut2Name] >= cut2)][self.parameter]*QtoE,
                            bins=binz,
                            range=(0,xhigh))[0]

    def get_histogram_2cuts_S_1_above_1_below(self,
                                              cut1Name,
                                              cut2Name,
                                              cut1,
                                              cut2,
                                              xhigh,
                                              binz) :
        
        return np.histogram(self.michelDF[(self.michelDF[cut1Name] <= cut1 ) & (self.michelDF[cut2Name] >= cut2)][self.parameter]*QtoE,
                            bins=binz,
                            range=(0,xhigh))[0]

    def get_histogram_2cuts_B_2_below(self,
                                      cut1Name,
                                      cut2Name,
                                      cut1,
                                      cut2,
                                      xhigh,
                                      binz) :
        
        return np.histogram(self.muonDF[(self.muonDF[cut1Name] <= cut1 ) & (self.muonDF[cut2Name] <= cut2)][self.parameter]*QtoE,
                            bins=binz,
                            range=(0,xhigh))[0]

    def get_histogram_2cuts_S_2_below(self,
                                      cut1Name,
                                      cut2Name,
                                      cut1,
                                      cut2,
                                      xhigh,
                                      binz) :
        
        return np.histogram(self.michelDF[(self.michelDF[cut1Name] <= cut1 ) & (self.michelDF[cut2Name] <= cut2)][self.parameter]*QtoE,
                            bins=binz,
                            range=(0,xhigh))[0]
