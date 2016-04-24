# Takes folder and move_threshold and print list of files in that folder which
# have 2d basecall and max_move <= move_threshold.
import h5py
import sys
import os
import argparse

description = """ Takes folder and move_threshold and print list of files in
that folder which have 2d basecall and max_move <= move_threshold."""

parser = argparse.ArgumentParser(description=description, epilog='')
parser.add_argument('--reads_folder', action='store',
                    help='Folder containing all the reads')
parser.add_argument('--move_threshold', action='store', type=int,
                    help="""Move move_threshold used for filtering reads.
                    We use only reads with max_move <= move_threshold.""")
parser.add_argument('--strand', action='store', help='Strand.')
args = parser.parse_args()

reads_folder = args.reads_folder
move_threshold = args.move_threshold
strand = args.strand

# End of parsing arguments.

if not reads_folder.endswith('/'):
    reads_folder += '/'

basecall_2d_node = "Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"
events_node = "Analyses/Basecall_2D_000/BaseCalled_%s/Events" % strand
for filename in os.listdir(reads_folder):
    if filename.endswith(".fast5"):
        try:
            fast5_file = h5py.File(reads_folder + filename, 'r')
            if basecall_2d_node in fast5_file:
                max_move = max(
                    map(lambda e: e["move"], fast5_file[events_node]))
                if max_move <= move_threshold:
                    print(reads_folder + filename)
        except Exception:
            pass
