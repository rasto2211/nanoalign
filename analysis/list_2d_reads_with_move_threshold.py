# Takes folder and move_threshold and print list of files in that folder which
# have 2d basecall and max_move <= move_threshold.
import h5py
import sys
import os

folder = sys.argv[1]
move_threshold = int(sys.argv[2])
strand = sys.argv[3]

if not folder.endswith('/'):
    folder += '/'

basecall_2d_node = "Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"
events_node = "Analyses/Basecall_2D_000/BaseCalled_%s/Events" %strand
for filename in os.listdir(folder):
    if filename.endswith(".fast5"):
        try:
            fast5_file = h5py.File(folder + filename, 'r')
            if basecall_2d_node in fast5_file:
                max_move = max(
                    map(lambda e: e["move"], fast5_file[events_node]))
                if max_move <= move_threshold:
                    print(folder + filename)
        except Exception:
            pass
