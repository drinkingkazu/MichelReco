import sys
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from ROOT import gSystem,TMath
from ROOT import michel
from ROOT import larlite as fmwk
from ROOT import larutil

##########################
def getRMS(vec,start,end):
    
    if (end <= start):
        print 'FAIL!'
        return 0
    if (start < 0):
        print 'FAIL!'
        return 0
    if (end >= len(vec)):
        print 'FAIL!'
        return 0

    arr = np.array(vec[start:end])
    base = np.mean(arr)
    rms  = np.std(arr)

    return base,rms

################
def getMax(vec):

    idx = np.argmax(vec)
    return idx,vec[idx]


##############################
def getDistToPoint(dS,peak,d):
    lpeak = dS[peak]
    npts = len(dS)-1
    for x in xrange(peak):
        l = lpeak-dS[peak-x]
        if (l > d):
            return peak-x
    # if we could not find a point
    # a distance d away from
    # the enf of dS
    print 'CANNOT FIND SUCH POINT!'
    return -1


###########################################
def getMIPIndices(truncQ,MIPmedian,MIPrms,cutoff):

    indices = []
    for n in xrange(cutoff):
        
        adc = truncQ[n]
        
        if ( (adc < MIPmedian+MIPrms) and (adc > MIPmedian-MIPrms) ):
            indices.append(n)

    return np.array(indices)

######################################
def getSubsetVectors(dS,truncQ,points):
    subS = []
    subQ = []
    for i in points:
        subS.append(dS[i])
        subQ.append(truncQ[i])

    return np.array(subS),np.array(subQ)


#################################
## get rms w.r.t. a linear slope
def getCorrectedRMS(f,MIPdS,MIPdQ):

    # line serves as the baseline
    rms  = 0
    for x in xrange(len(MIPdS)):
        rms += (MIPdQ[x]-f(MIPdS[x])) * (MIPdQ[x]-f(MIPdS[x]))
        rms = np.sqrt(rms)

    return rms

##########################
# get area from MIP end to bragg peak accounting for slope
def getArea(f,dS,truncQ,MIPend,maxIdx):

    area = 0
    for x in xrange(MIPend,maxIdx+1):
        area += truncQ[x]-f(dS[x])

    return area

# Create ana_processor instance
my_proc = fmwk.ana_processor()

my_proc.enable_filter(True)
the_filter = fmwk.StoppingFilter()
the_filter.Filter(False)
the_filter.Stopping(False)

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify analysis output root file name
my_proc.set_ana_output_file("out.root");

plot = False
_verb = False

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])


stopping = fmwk.TagStoppingMuon()
stopping.SetClusterProducer('cccluster')
stopping.SetAlgo(michel.TruncatedQBoundary())

my_proc.add_process(the_filter)
my_proc.add_process(stopping)

print
print  "Finished configuring ana_processor. Start event loop!"
print

braggsig = []
muE = []
muEstop = []
muEesc  = []
braggDist = []
braggamp = []
braggarea = []
braggareaSTOP = []
braggareaESC  = []
bragglen = []

if (plot):
    fig, ax = plt.subplots(figsize=(10,8))
    plt.ion()

while (my_proc.process_event()):

    if (the_filter.Passes() == False):
        continue


    truncQ = np.array(stopping.GetTruncQ())
    dS     = np.array(stopping.GetdS())

    if (len(truncQ)):

        # if < 30 cm -> continue
        if (dS[-1] < 30):
            continue

        maxIdx, maxVal = getMax(truncQ)
        
        # figure out which side the maxIdx (i.e. position of bragg peak is at)
        if (maxIdx < len(truncQ)/2.):
            dS = dS-dS[-1]
            truncQ = truncQ[::-1]
            maxIdx = len(dS)-maxIdx-1

        # focus on MIP portion for now
        # skip the last 15 cm till bragg peak
        # get the RMS and median value
        MIPend = getDistToPoint(dS,maxIdx,15)
        MIPmedian = np.median(np.array(truncQ[:MIPend]))
        MIPrms    = np.std(np.array(truncQ[:MIPend]))
        if (_verb):
            print 'MIP median: %.02f\tMIP rms: %.02f'%(MIPmedian,MIPrms)
        # select hits within 1 RMS of the MIPmedian
        # get list of indices for vector truncQ that match this req
        MIPindices = getMIPIndices(truncQ,MIPmedian,MIPrms,MIPend)
        if (len(MIPindices) == 0):
            continue
        MIPdS,MIPdQ = getSubsetVectors(dS,truncQ,MIPindices)
        p = np.polyfit(MIPdS,MIPdQ,1)
        f = np.poly1d(p)
        # now that we have a function to model the slope
        # calculate RMS of MIP points w.r.t this function
        corrRMS = getCorrectedRMS(f,MIPdS,MIPdQ)
        if (_verb):
            print 'corrected RMS: %.02f'%corrRMS

        # now calculate the "significance" of the Bragg peak
        # w.r.t. baseline at that dS value
        braggPeakExpected = f(dS[maxIdx])
        braggAmp = maxVal-braggPeakExpected
        braggSignificance = braggAmp/corrRMS
        if (_verb):
            print 'Bragg dS           : %.02f'%dS[maxIdx]
            print 'Expected Bragg Peak: %.02f'%braggPeakExpected
            print 'Bragg Peak         : %.02f'%maxVal
            print 'Bragg Ampltidue    : %.02f'%braggAmp
            print 'Bragg Significance : %.02f'%braggSignificance
        braggamp.append(braggAmp)
        #braggsig.append(braggSignificance)
        braggsig.append(MIPrms)

        # distance of bragg peak to end of track
        dist = np.abs(dS[maxIdx]-dS[-1])
        # distance of bragg peak to end of MIP
        braggLen = np.abs(dS[maxIdx]-dS[MIPend])
        if (_verb):
            print 'Bragg Dist to End  : %.02f'%dist
            print 'Bragg-MIP end dist : %.02f'%braggLen
        braggDist.append(dist)
        bragglen.append(braggLen)

        # area up to bragg peak
        braggArea = getArea(f,dS,truncQ,MIPend,maxIdx)
        if (_verb):
            print 'Area is            : %.02f'%braggArea
        braggarea.append(braggArea)

        E = the_filter.GetEndEnergy()
        muE.append(E)       
        if (E < 10):
            braggareaSTOP.append(braggArea)
            muEstop.append(E)
        elif (E > 10):
            braggareaESC.append(braggArea)
            muEesc.append(E)
        if (_verb):
            print "End energy is: %.02f"%E


        if (plot):

            ax.cla()
            fig.gca()            
            ax.plot(dS,truncQ,'bo--')
            plt.scatter(dS[maxIdx],maxVal,color='r',s=300)
            plt.plot(MIPdS,MIPdQ,'mo',markersize=10,alpha=0.5)
            plt.scatter(dS[MIPend],truncQ[MIPend],color='g',s=300)
            plt.plot(MIPdS,f(MIPdS),'k-',lw=3)
            plt.grid()
            plt.ylim([0,250])
            plt.show()
            
            usrinput = raw_input("Hit Enter: next evt  || int: go to event number ||  q: exit viewer\n")                        
            if ( usrinput == "q" ):
                sys.exit(0)
            else:
                continue

#plt.scatter(braggsig,muE)
#plt.scatter(muE,braggDist)
#plt.scatter(muE,braggamp)
#plt.scatter(muE,braggarea)
plt.hist2d(muEesc,braggareaESC,bins=(np.linspace(0,500,40),np.linspace(-1300,3000,40)))
#plt.xlim([-20,100])
#plt.ylim([0,1000])
plt.grid()
plt.show()

plt.hist(braggareaSTOP,np.linspace(-1300,3300,30),color='r',label='Stopping',alpha=0.5,normed=True)
plt.hist(braggareaESC,np.linspace(-1300,3300,30),color='b',label='Escaping',alpha=0.5,normed=True)
plt.grid()
plt.legend()
plt.show()

sys.exit(0)

