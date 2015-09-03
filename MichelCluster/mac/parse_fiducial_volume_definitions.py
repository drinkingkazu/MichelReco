import os
import ROOT

textfile = os.environ['LARLITE_USERDEVDIR'] + '/MichelReco/MichelCluster/mac/fiducial_volume_definitions.txt'

if not os.path.isfile(textfile):
    print "You're using parse_fiducial_volume_definitions.py and it's looking for a file that doesn't exist!"
    print "The file it wants is",textfile
    print "Quitting!"
    quit()

def list_wires_times_to_exclude():
    parsing_wires, parsing_times = False, False
    wires_to_exclude_min, times_to_exclude_min = ROOT.std.vector("double")(), ROOT.std.vector("double")()
    wires_to_exclude_max, times_to_exclude_max = ROOT.std.vector("double")(), ROOT.std.vector("double")()
    for line in open(textfile):
        if line[:13] == 'exclude wires': 
            parsing_wires = True
            parsing_times = False
            continue
        elif line[:13] == 'exclude times':
            parsing_wires = False
            parsing_times = True
            continue
        #skip any lines that don't start with an integer
        try: 
            int(line[0])
        except: pass
    
        #skip empty lines
        if len(line.strip('\n').split('-')) < 2: continue
        
        region_to_exclude = tuple(line.strip('\n').split('-'))
        region_to_exclude_min = float(region_to_exclude[0])
        region_to_exclude_max = float(region_to_exclude[1])

        if parsing_wires and not parsing_times: 
            print 'Removing Wire region[%i,%i]'%(int(region_to_exclude_min),int(region_to_exclude_max))
            wires_to_exclude_min.push_back(region_to_exclude_min)
            wires_to_exclude_max.push_back(region_to_exclude_max)
        elif not parsing_wires and parsing_times:
            print 'Removing Time region[%i,%i]'%(int(region_to_exclude_min),int(region_to_exclude_max))
            times_to_exclude_min.push_back(region_to_exclude_min)
            times_to_exclude_max.push_back(region_to_exclude_max)
        else:
            print 'wtf kaleko screwed something up in the script that parses fiducial volume info'

    return wires_to_exclude_min, wires_to_exclude_max, times_to_exclude_min, times_to_exclude_max
