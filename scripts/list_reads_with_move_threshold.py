# Reads filtered_by_stats.csv from stdin and filter them by move and output
# final files to stdout.
import h5py
import sys

move_threshold = int(sys.argv[1])

events_node = "Analyses/Basecall_2D_000/BaseCalled_template/Events"
header = input()
for line in sys.stdin:
    try:
        path = line.split(',')[0].strip('"')[:-len(".sam")] + ".fast5"
        fast5_file = h5py.File(path, 'r')
        max_move = max(map(lambda e: e["move"], fast5_file[events_node]))
        if max_move <= move_threshold:
            print(path)
    except Exception:
        pass
