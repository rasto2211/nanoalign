# Takes SAM file with alignment of the read with reference sequence and returns
# the part of ref. seq. that was aligned to the read.
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

filename = sys.argv[1]
file = pysam.Samfile(filename)

for read in file:
    cigartuples = read.cigartuples
    if cigartuples is None:
        print("NA")
        continue

    cigar_counts = cigar_profile(cigartuples)
    # insertions + deletions + matches + mismatches
    MID = cigar_counts["INS"] + cigar_counts["DEL"] + cigar_counts["MATCH"]

    edit_distance = read.get_tag("NM")
    identity = (1 - edit_distance / MID) * 100
    print(">%s Identity: %f" % (filename, identity))
    print(read.get_reference_sequence().upper())
