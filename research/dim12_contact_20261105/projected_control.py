"""A separate 841 discovery control for the prior projected-Adam explorer.

This explicitly differs from the source's unprojected float64 Adam search.
It uses the same four families as the prior 842 exploratory campaign.
"""
from common import *
sys.path.insert(0,str(AUDIT))
import riesz_search as explorer

original_initialize=explorer.initialize
def initialize(batch):
    x,metadata=original_initialize(batch);x=x[:,:841].copy()
    for row,y in zip(metadata,x):
        row['initial']=metrics(y);row['control_count']=841
    return x,metadata

if __name__=='__main__':
    explorer.initialize=initialize
    explorer.main()
