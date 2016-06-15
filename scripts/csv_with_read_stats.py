# Produces CSV file with the stats taken from all the SAM files in the
# directory given as a second argument. The first argument is the number of jobs
# in multiprocessing pool (optional). If the second argument is not present 
# then the list of files is read from stdin.

import sys
import sam_utils
import os

from multiprocessing import Pool


def get_read_stats(file_path):
    alignment = sam_utils.get_best_alignment(file_path)
    if alignment is None:
        return None

    read_name = file_path[:-len(".sam")]

    fasta_file = open(read_name + ".fasta")
    fasta_header = fasta_file.readline().strip()
    whole_read = fasta_file.readline().strip()

    cigar_counts = sam_utils.get_cigar_counts(alignment)
    edit_distance = sam_utils.get_edit_distance_from(alignment)

    return (read_name, edit_distance,
            cigar_counts["INS"], cigar_counts["DEL"],
            cigar_counts["MATCH"], cigar_counts["SOFT_CLIP"],
            cigar_counts["HARD_CLIP"], cigar_counts["SKIP"],
            cigar_counts["PAD"], alignment.query_alignment_start,
            alignment.query_alignment_end, alignment.query_length,
            len(whole_read), alignment.reference_start,
            alignment.reference_end, alignment.reference_length)

jobs = 1
if len(sys.argv) == 2:
    jobs = int(sys.argv[1])

file_paths = []
if len(sys.argv) == 3:
    path_to_dir = sys.argv[2]
    file_paths = [path_to_dir + "/" +
                  file for file in os.listdir(path_to_dir) if file.endswith(".sam")]
else:
    file_paths = [path.strip() for path in sys.stdin]

pool = Pool(jobs)
chunk = len(file_paths) // jobs
stats = pool.map(get_read_stats, file_paths, chunk)
pool.close()

csv_header = str("path, edit_dist, I, D, M, S_CLIP, H_CLIP, SKIP, PAD, "
                 "query_start, query_end, query_length, read_len, ref_start, ref_end,"
                 " ref_length")
print(csv_header)
for row in stats:
    if row is None:
        cols = len(csv_header.split(','))
        print(",".join(["NA"] * cols))
        continue

    print(",".join(map(str, list(row))))
