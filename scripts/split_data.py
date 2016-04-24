import h5py
import sys
import os
import random
import argparse
import datetime

description = """
Splits data into training and testing set based on number of events.
Results are written into training_set.list and testing_set.list.
"""
parser = argparse.ArgumentParser(description=description, epilog='')
parser.add_argument('--file_with_list', action='store',
                    help='File containing list to be used.')
parser.add_argument('--strand', action='store', help='Strand')
parser.add_argument(
    '--training_set_percent', action='store', type=float, default=0.7,
    help='How much data will be used fo training? (% of all events)')
args = parser.parse_args()

########################################################################

file_with_list = args.file_with_list
strand = args.strand
training_set_percent = args.training_set_percent

file_paths = [path.strip() for path in open(file_with_list)]

random.seed(datetime.datetime.now())
random.shuffle(file_paths)

events_node = "Analyses/Basecall_2D_000/BaseCalled_%s/Events" % strand
path_filesize = [(path, len(h5py.File(path, 'r')[events_node]))
                 for path in file_paths]

total_size = sum(f[1] for f in path_filesize)

training_set_size = 0
training_set = open("training_set.list", 'w')
testing_set = open("testing_set.list", 'w')
for path, size in path_filesize:
    if (training_set_size + size) / total_size <= training_set_percent:
        training_set.write("%s\n" % path)
        training_set_size += size
    else:
        testing_set.write("%s\n" % path)

print("Total events: %d" %(total_size))
print("Training set events: %d" %training_set_size)
