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
    wires_to_exclude, times_to_exclude = ROOT.std.vector(ROOT.std.pair("double","double"))(), ROOT.std.vector(ROOT.std.pair("double","double"))()
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
        region_to_exclude = ROOT.std.pair("double","double")(float(region_to_exclude[0]),float(region_to_exclude[1]))

        if parsing_wires and not parsing_times: wires_to_exclude.push_back(region_to_exclude)
        elif not parsing_wires and parsing_times: times_to_exclude.push_back(region_to_exclude)
        else: print 'wtf kaleko screwed something up in the script that parses fiducial volume info'
        
    return wires_to_exclude, times_to_exclude
