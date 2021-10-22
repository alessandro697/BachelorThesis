
#!/usr/bin/env python3
import sys
import re
import os
import os.path
import numpy as np # to install this lib in deb-type linux, run as root: apt-get install python3-numpy
import xyz_utils
from collections import namedtuple

def convert_1Dfiles2xyz(filenames,printout=False):
    """read output file from 1D FK, in format time x1 x2.... v1 v2...
    and generate output in xyz format"""
    if len(filenames)<1:
        print("error: need at least 1 file!!!", file=sys.stderr)
        return None
    
    allfr=[]
    for filen in filenames:
        if filen=="-":
            f=sys.stdin
        else:
            f = open(filen, 'r')
        for line in f:
            line=re.sub("^\s+","",line).rstrip()
            nums=re.split('\s+',line)
            time=nums[0]
            comment="# time= "+time+" from file "+filen
            natom=(len(nums)-1)//2
            atoms=['H']*natom
            coords = np.zeros([natom, 3], dtype="float")
            for i in range(natom):
                coords[i,0]=float(nums[i+1])
            fr=namedtuple("xyzdata", ["comment", "atoms", "coords"])(comment, atoms, coords)
            if printout:
                xyz_utils.xyz_write_one_frame(sys.stdout,fr)
            else:
                allfr.append(fr) # memory-dangerous! cumulating all frames in memory...
    return allfr
        

# the following function is only exectuted when this code is run as a script, and its purposes is to parse
# the command line and to generate a meaningful parsed args list to the actual function doing the job:
if __name__ == "__main__":
    import sys
    import argparse
    commandname=sys.argv[0]
# default values:
    filenames = []

    desc="""convert 1D output from fk.py to a standard xyz format
    INPUT: xyz files
    OUTPUT: a regular xyz file"""

    epil="""			v1.0 by N. Manini, 02.05.2020"""

    ##  Argument Parser definition: this is just an example...
    parser = argparse.ArgumentParser( formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                    , description=desc, epilog=epil)

    parser.add_argument( 'filenames', nargs='*', default=['-'],
                         help='Files to be processed. If not given, stdin is used')
   
    ## End arg parser definition
    args=parser.parse_args(sys.argv[1:])
    d = vars(args)	# adding prog in args, for unknown reasons it's not there...
    d['prog']=parser.prog
#   here the actual function doing the job is called:
    a=convert_1Dfiles2xyz(args.filenames,True)
#   final printout:
#	xyz_utils.xyz_write_entire_file(sys.stdout, a)

