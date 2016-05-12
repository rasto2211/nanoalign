# Reads filtered_by_stats.csv from stdin and filter them by move and output
# final files to stdout.
import h5py
import sys

from multiprocessing import Pool


def max_move(path):
    try:
        events_node = "Analyses/Basecall_1D_000/BaseCalled_template/Events"
        fast5_file = h5py.File(path, 'r')
        return max(map(lambda e: e["move"], fast5_file[events_node]))
    except:
        events_node = "Analyses/Basecall_2D_000/BaseCalled_template/Events"
        fast5_file = h5py.File(path, 'r')
        return max(map(lambda e: e["move"], fast5_file[events_node]))

jobs = int(sys.argv[1])
move_threshold = int(sys.argv[2])

header = input()
paths = [line.split(',')[0].strip('"') + ".fast5" for line in sys.stdin]

pool = Pool(jobs)
chunk = len(paths) // jobs
max_moves = pool.map(max_move, paths, chunk)

for move, path in zip(max_moves,paths):
    print("%d %s" %(move,path))

pool.close()
