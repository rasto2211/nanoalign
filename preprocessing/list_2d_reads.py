# This script lists all paths to 2D reads that are in the folder that is given
# in arguments. It can be used as input for training HMM.
import h5py
import sys
import os

folder = sys.argv[1]

if not folder.endswith('/'):
    folder += '/'

basecall_2d = "Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"
for filename in os.listdir(folder):
    if filename.endswith(".fast5"):
        try:
            fast5_file = h5py.File(folder + filename, 'r')
            if basecall_2d in fast5_file:
                print(filename)
        except Exception:
            pass
