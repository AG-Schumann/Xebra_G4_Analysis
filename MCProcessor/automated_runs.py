
#!/usr/bin/python

import sys
import os
from os import listdir
from os.path import isfile, join

'''
Automated clustering of all ROOT files in input folder using DBScan.
'''

## Create list of all files in input folder
PATH = 'input/'
files = [f for f in listdir(PATH) if isfile(join(PATH, f))]

## Iteratively loop over all files in aforementioned list to be clustered
for FILENAME in files:
    
    ## Check if ROOT file
    if FILENAME.find('root') == -1:
        continue
        
    ## Remove '.root' extension
    FILENAME = FILENAME[:-5]
    
    ## Print feedback in terminal
    sys.stdout.write("\nStarted processing of "+PATH+FILENAME+".root\n\n")
            
    # Execute clustering
    os.system("./automated_runs.sh %s %s"%(PATH, FILENAME))
            
sys.stdout.write("\nDone.\n")

