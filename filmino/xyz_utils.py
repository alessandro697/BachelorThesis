#!/usr/bin/env python3

# loosely inspired by
# https://github.com/pele-python/pele/blob/master/pele/utils/xyz.py

import re
from collections import namedtuple
import numpy as np # to install this lib in deb-type linux, run as root: apt-get install python3-numpy

def xyz_read_one_frame(f):
    """read and return one frame from a xyz file
    f should be a filehandle to a previously opened xyz file"""
    debug=False
    natom=0
    count=-2
    comment=None
    atoms=None
    coords=None
    for line in f:
        line=re.sub("^\s+","",line).rstrip()
        nums=re.split('\s+',line)
        if debug:
            print(" debug ", len(nums), " nums[0]: ",nums[0]," coords: ", coords)
        if len(nums)>=4 and count>=0:
            atoms.append(nums[0])
            coords[count,:]=list(map(float, nums[1:4]))
            count+=1
        elif count == -1:
            count=0
            comment=line
            if debug:
                print ("debug: comment=",line)
        elif count==-2 and len(nums)==1:
            count=-1
            natom=int(nums[0])
            atoms=[]
            coords = np.zeros([natom, 3], dtype="float")
        else:
            print(commandname, "ERROR, invalid line in xyz file:\n",line, file=sys.stderr)
#            break
        if debug:
            print ("debug: count=",count, natom)
        if count>=natom: # have read enough xyz data: frame is complete, can terminate here
            break
    return namedtuple("xyzdata", ["comment", "atoms", "coords"])(comment, atoms, coords)


def xyz_read_entire_file(f):
    """read and return an entire xyz file
    f should be a filehandle to a previously opened xyz file"""
    allf=[]
    while True: # keeps looping until the end of the file
        snap=xyz_read_one_frame(f)
        if snap.comment==None:
            break
        else:
            allf.append(snap) # this is potentially quite dangerous, "allf" may suck a LOT of memory!
    return allf


def xyz_write_one_frame(fout, xyzdata):
    """append a single frame to a xyz file
    fout should be a filehandle to a writable file (e.g. sys.stdout is fine)"""
                                        # xyzdata is a namedtuple("xyzdata", ["comment", "atoms", "coords"])(comment, atoms, coords)
    fout.write("%d\n%s\n" % (len(xyzdata.atoms), xyzdata.comment))
#    print(len(xyzdata.atoms),file=fout)
#    print(xyzdata.comment,file=fout)
    for i in range(len(xyzdata.atoms)):
        fout.write("%s %g %g %g\n" % (xyzdata.atoms[i], xyzdata.coords[i,0], xyzdata.coords[i,1], xyzdata.coords[i,2]))

def xyz_write_entire_file(fout, m_xyz):
    """append multiple frames to a xyz file
    fout should be a filehandle to a writable file (e.g. sys.stdout is fine)"""

                                        # m_xyz is a list of namedtuple("xyzdata", ["comment", "atoms", "coords"])(comment, atoms, coords)
    for snap in m_xyz:
        xyz_write_one_frame(fout,snap)


# the following function is only exectuted when this code is run as a script, and its purposes is to parse
# the command line and to generate a meaningful parsed args list to the actual function doing the job:
if __name__ == "__main__":
    import sys
    import argparse
    commandname=sys.argv[0]

    desc="""read and write xyz files.
    This function serves just for testing the functions provided by this library
    """

    epil="""		v. 1.2	by Nicola Manini, 23/04/2020"""

    parser = argparse.ArgumentParser( formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                    , description=desc, epilog=epil)

    parser.add_argument( 'filenames', nargs='*', default=['-'],
                         help='xyz files to be read. If none is given, stdin is used')

#    parser.add_argument( '-d', '--debug', action='store_true',
#                         dest='debug', 
#                         help='activate debug mode -- WARNING! output is affected/spoiled!' )

    args=parser.parse_args(sys.argv[1:])
    d = vars(args)	# adding prog in args, for unknown reasons it's not there...
    d['prog']=parser.prog

    for filen in args.filenames:
        print("#now readinging",filen,"one frame at a time")
        count=0
        if filen=="-":
            f=sys.stdin
        else:
            f = open(filen, 'r')

        nframe=0
#   here the library functions are tested:
        while True:
            snap=xyz_read_one_frame(f)
            if snap.comment==None:
                break
            nframe+=1
            print("#this is frame",nframe, "of file",filen)
#   the results are printed out:
            print(len(snap.atoms))
            print(snap.comment)
            for i in range(len(snap.atoms)):
#                print(snap.atoms[i],snap.coords[i])
                print(snap.atoms[i]," ".join([str(snap.coords[i,j]) for j in range(3)]))
        if filen!="-":
            f.close()
        print("Read and processed correctly",nframe,"frames\n\n")

        print("test reading the entire file",filen,"in one shot:")
        if filen=="-":
            f=sys.stdin
        else:
            f = open(filen, 'r')
        allf=xyz_read_entire_file(f)
        print(allf)
        if filen!="-":
            f.close()
        print("Read and processed correctly",len(allf),"frames, that formally occupy",sys.getsizeof(allf),"bytes\n\n")

        print("Writing in sys.stdout the first frame of those previously loaded")
        xyz_write_one_frame(sys.stdout, allf[0])

        filen="/tmp/test1.xyz";
        print("Writing in ",filen,"first max 2 frames of those previously loaded")
        
        f = open(filen, 'w')
        xyz_write_entire_file(f, allf[:2])
        f.close()

        print("Done writing in ",filen,"\n\n")
