# This script prints statistics from SAM file which was produced after alignment
# of two sequences by BWA-MEM.
import sys
import sam_utils

print("edit_distance\tMID\tquery_start\tquery_end\tref_start\tref_end")
for alignment in sam_utils.get_all_alignments_from(sys.argv[1]):
    print("\t".join(map(str, list(alignment))))
