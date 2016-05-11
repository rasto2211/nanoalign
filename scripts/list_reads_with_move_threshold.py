# Reads filtered_by_stats.csv from stdin and filter them by move and output
# final files to stdout.
import h5py
import sys
import functools

from multiprocessing import Pool


def exceeds_move_threshold(path, move_threshold):
    events_node = "Analyses/Basecall_2D_000/BaseCalled_template/Events"
    fast5_file = h5py.File(path, 'r')
    max_move = max(map(lambda e: e["move"], fast5_file[events_node]))
    if max_move > move_threshold:
        return True

    return False

jobs = int(sys.argv[1])
move_threshold = int(sys.argv[2])

header = input()
paths = [line.split(',')[0].strip('"') + ".fast5" for line in sys.stdin]

pool = Pool(jobs)
chunk = len(paths) // jobs
exceeds_threshold = pool.map(
    functools.partial(exceeds_move_threshold, move_threshold=move_threshold),
    paths, chunk)
pool.close()

for exceeds, path in zip(exceeds_threshold, paths):
    if not exceeds:
        print(path)
