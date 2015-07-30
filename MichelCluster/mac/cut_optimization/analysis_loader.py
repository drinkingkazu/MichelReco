import ROOT
import root_numpy
import pandas as pd


#will be signal or background :)

def get_data_frame(tfile,ttree):
    rec = root_numpy.root2array(tfile,ttree)
    return pd.DataFrame(rec)
