import ROOT
from ROOT import michel
import parse_fiducial_volume_definitions as fidparser

def PrepareMichelAlgo():

    algoList = []


    ########################################
    # remove clusters w/ less than some #
    # of hits
    minclushits = michel.CutOnTotNumHits()
    minclushits.SetMinTotHits(25)
    #!algoList.append(minclushits)

    #########################################
    # calculate various cluster parameters...
    ctrunk = michel.CalcTruncated()
    #ctrunk.SetCovarianceWindowSize(int s)      
    #ctrunk.SetTruncatedQWindowSize(int s)      
    #ctrunk.SetPAbove(double p)                 
    #ctrunk.SetMinWindowSize(int w)             
    #ctrunk.SetEdgeEffectFix(int e)             
    algoList.append(ctrunk)

    ########################################
    # Attach algorithm for boundary finding
    boundaryalgo = michel.BoundaryFromTQMaxQ()
    boundaryalgo.SetMaxDistancesTruncatedQMaxQ(15)
    algoList.append(boundaryalgo)

    ##############################################
    # Attach algo for charge spectrum requirements
    closepeaks = michel.RequireCloseTruncatedPeaks()
    closepeaks.SetMaxDistanceTruncatedPeaks(5)
    #algoList.append(closepeaks)

    covdip = michel.RequireCovarianceDip()
    covdip.SetCovarianceDipCutoff(0.9)
    #algoList.append(covdip)

    slopeflip = michel.RequireSlopeSignFlip()
    #algoList.append(slopeflip)

    lowcovbound = michel.RequireBoundaryInLowCov()
    lowcovbound.SetMaxCovarianceAtStart(0.8)
    algoList.append(lowcovbound)

    #############################################
    # Attach algorithm for finding michel cluster
    findMichel = michel.ForwardMichelID()
    findMichel.SetMaxMichelHits(0)
    algoList.append(findMichel)

    #########################
    # MID finding algorithms
    midalgo = michel.DecideIfStoppingMuon()
    midalgo.SetChiMin     ( 0.9 )
    midalgo.SetFracMinHits( 0.7 )
    midalgo.SetHitRadius  ( 30  )
    midalgo.SetMaxDist    ( 3.0 )
    midalgo.SetMinBadHits ( 10  )
    algoList.append(midalgo)

    #########################
    # BraggArea filter algo
    braggalgo = michel.FindBraggPeak()
    braggalgo.SetMinBraggArea(1000.)
    algoList.append(braggalgo)

    ############################
    # MINIMUM LENGTH REQUIREMENT
    minlength = michel.CutOnMuonLength()
    minlength.SetMinMuonLength(10)
    algoList.append(minlength)

    ############################
    # MUON LINEARITY CUT
    minlinearity = michel.CutOnMuonLinearity()
    minlinearity.SetChiMin     ( 0.8 )
    minlinearity.SetFracMinHits( 0.5 )
    algoList.append(minlinearity)

    #########################################################
    # MID filter that removes michels close to wire gaps/edges
    fidvolfilter = michel.CutOnFiducialVolume()

    wires_to_exclude_min, \
        wires_to_exclude_max, \
        times_to_exclude_min, \
        times_to_exclude_max = fidparser.list_wires_times_to_exclude()
    
    fidvolfilter.SetExcludedWireRanges(wires_to_exclude_min,wires_to_exclude_max)
    fidvolfilter.SetExcludedTimeRanges(times_to_exclude_min,times_to_exclude_max)
    #algoList.append(fidvolfilter)

    ######################################
    # Attach algorithm to recluster michel
    supersonic = michel.SuperSonicClusterer()
    supersonic.SetMergeTillConverge(True)
    supersonic.SetUseHitRadius     (True)
    supersonic.SetMaxRadius( 15 )
    supersonic.SetHitRadius(  3 )
    #supersonic.SetVerbosity(michel.msg.kDEBUG)
    algoList.append(supersonic)

    #########################
    # cone-finding algorithm
    conefinder = michel.ConeHitFinder()
    conefinder.SetMaxRadius               ( 20 )
    conefinder.SetMaxPerpendicularDistance(  3 )
    algoList.append(conefinder)

    ##############################################
    # final mid algo cutting on num of michel hits
    michelhits = michel.CutOnMichelNumHits()
    michelhits.SetMinMichelHits (  5 )
    michelhits.SetMaxMichelHits ( 35 )
    algoList.append(michelhits)

    #########################################
    # remove weird horizontal tracks from PMT
    pmtremoved = michel.RemoveFakePMTSignals()
    pmtremoved.SetMaxErrorTime(0.1)
    algoList.append(pmtremoved)

    #############################################
    # require large angle between michel and muon
    largeangle = michel.RequireLargeAngle()
    largeangle.SetMinAngle(30.*3.14/180.)
    largeangle.SetMinStraightMichelHits(5)
    largeangle.SetMuonLengthUsed(10000)
    #largeangle.SetVerbosity(michel.msg.kDEBUG)
    algoList.append(largeangle)

    #############################################
    # remove bragg peak hits
    removeBragg = michel.RemoveBraggPeakHits()
    removeBragg.SetMaxRadius(1.)
    removeBragg.SetChargeFactor(2.)
    algoList.append(removeBragg)


    #####################################
    # cut on Michels with high avg. Q/hit
    cutonavgq = michel.CutOnMeanHitCharge()
    cutonavgq.SetMaxAvgQ(400.)
    algoList.append(cutonavgq)

    return algoList
