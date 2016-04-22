# This script prints statistics from SAM file which was produced after alignment
# of two sequences by BWA-MEM.
import sys
import pysam

from collections import Counter

code2str = {
    0: "MATCH",  # M - alignment match (can be a sequence match or mismatch)
    1: "INS",  # I
    2: "DEL",  # D
    3: "SKIP",  # N
    4: "SOFT_CLIP",  # S
    5: "HARD_CLIP",  # H
    6: "PAD",  # P
    7: "EQUAL",  # =
    8: "DIFF",  # X
}


def cigar_profile(cigar_tuples):
    cigar_prof = Counter()
    for cigar_tuple in cigar_tuples:
        cigar_operation = code2str[cigar_tuple[0]]
        cigar_prof[cigar_operation] += cigar_tuple[1]
    return cigar_prof

file = pysam.Samfile(sys.argv[1])

print("edit_distance\tMID\tquery_start\tquery_end\tref_start\tref_end")
for read in file:
    cigartuples = read.cigartuples
    if cigartuples is None:
        print("NA")
        continue
    cigar_counts = cigar_profile(cigartuples)
    # insertions + deletions + matches + mismatches
    MID = cigar_counts["INS"] + cigar_counts["DEL"] + cigar_counts["MATCH"]

    query_start = read.query_alignment_start
    query_end = read.query_alignment_end
    query_len = read.query_alignment_length
    ref_start = read.reference_start
    ref_end = read.reference_end
    ref_len = read.reference_length
    edit_distance = read.get_tag("NM")
    print("%d\t%d\t%d\t%d\t%d\t%d" %
          (edit_distance, MID, query_start, query_end, ref_start, ref_end))

    # print("query start %d end %d len %d" %(query_start,query_end,query_len))
    # print("ref start %d end %d len %d" %(ref_start,ref_end,ref_len))
    # print("edit distance  %d" %edit_distance)
    # print("inferred_length %d" %read.infer_query_length())
    # print("length of overlap %d" %read.get_overlap(start,end))
    # identity = (1 - edit_distance/MID)*100
    # print("identity %f %%" %identity)
