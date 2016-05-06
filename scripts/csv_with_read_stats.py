# Produces CSV file with the stats taken from all the SAM files in the 
# directory given as a second argument. The first argument is number of 
# jobs.
import sys
import sam_utils
import os

from multiprocessing import Pool


def get_read_stats(file_path):
    alignment, best_identity, percent_clipped = sam_utils.get_best_alignment(
        file_path)

    if alignment is None:
        return None

    cigar_counts = sam_utils.get_cigar_counts(alignment)

    return (file_path, best_identity, percent_clipped,
            cigar_counts["INS"], cigar_counts["DEL"], cigar_counts["MATCH"],
            cigar_counts["SOFT_CLIP"], cigar_counts["HARD_CLIP"],
            alignment.query_length,
            alignment.reference_length)

jobs = int(sys.argv[1])
path_to_dir = sys.argv[2]

pool = Pool(jobs)
file_paths = [path_to_dir + "/" +
              file for file in os.listdir(path_to_dir) if file.endswith(".sam")]
chunk = len(file_paths) // jobs
stats = pool.map(get_read_stats, file_paths, chunk)
pool.close()

print("path, id, clip, ins, del, match, "
"soft_clip, hard_clip, query_length, ref_length")
for row in stats:
    if row is None:
        continue

    print(",".join(map(str, list(row))))
