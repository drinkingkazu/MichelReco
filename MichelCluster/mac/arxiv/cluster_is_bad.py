import numpy as np

def isBad(all_hit_w,all_hit_t,clus_w,clus_t,clus_chi,clus_slope,boundary,forward,start_w,start_t):

    slope = 0
    count = 0
    w_avg = 0
    t_avg = 0

    print 'boundary is: ',boundary

    if (forward):
        for n in xrange(boundary):
            if (clus_chi[n] > 0.9):
                slope += clus_slope[n]
                count += 1
                w_avg += clus_w[n]
                t_avg += clus_t[n]
    else:
        for n in xrange(boundary,len(clus_slope)):
            if (clus_chi[n] > 0.9):
                slope += clus_slope[n]
                count += 1
                w_avg += clus_w[n]
                t_avg += clus_t[n]

    # require that a minimum fraction of the hits in muon be used
    frac = float(count)/float(len(clus_slope))
    print 'fraction of hits used to get slope: %.02f'%frac

    if (frac < 0.7):
        print 'ignoring because not enough strait points...'
        return [],[]

    slope /= count
    w_avg /= count
    t_avg /= count

    # get intercept (anchor) for line by using average point
    b = t_avg - slope * w_avg

    print 'Slope    : %.02f'%slope
    print 'Intercept: %.02f'%b
    print 'avg point: [%.02f,%.02f]'%(w_avg,t_avg)

    print 'Michel start w: %.02f'%start_w
    print 'Michel start t: %.02f'%start_t

    bad_w = []
    bad_t = []
    
    # identify bad hits from muon
    nbad = 0
    for h in xrange(len(all_hit_w)):

        hw = all_hit_w[h]
        ht = all_hit_t[h]

        # distance from michel start:
        dd = (start_w-hw) * (start_w-hw) + (start_t-ht) * (start_t-ht)
        
        if (dd < 100):

            print 'hit close is [%.02f,%.02f]'%(hw,ht)
            
            # check that this hit is not in cluster vector
            use = True
            for i in xrange(len(clus_w)):
                if ( (hw == clus_w[i]) and (ht == clus_t[i]) ):
                    use = False
                    break

            if (use == True):
                # get projection onto slope line
                proj = np.abs( ht - slope * hw - b ) / np.abs(slope)
                print 'proj: %.02f'%proj
                if (proj < 0) : proj *= -1
                if (proj < 3) :
                    nbad += 1
                    bad_w.append(hw)
                    bad_t.append(ht)


    print 'Bad points: %i'%nbad

    if (nbad > 10):
        return bad_w,bad_t

    return [],[]

