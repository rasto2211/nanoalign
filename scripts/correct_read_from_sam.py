# Takes SAM file with alignment of the read with reference sequence and returns
# the part of ref. seq. that was aligned to the read.
import sys
import pysam
import sam_utils

from collections import namedtuple

filename = sys.argv[1]
file = pysam.Samfile(filename)

Alignment = namedtuple("Alignment", "indentity ref filename")
best_alignment = None
for alignment in file:
    cigartuples = alignment.cigartuples
    if cigartuples is None:
        continue

    cigar_counts = sam_utils.cigar_profile(cigartuples)
    # insertions + deletions + matches + mismatches
    MID = sam_utils.get_MID_from(cigar_counts)

    edit_distance = sam_utils.get_edit_distance_from(alignment)
    identity = (1 - edit_distance / MID) * 100

    if best_alignment is None or identity > best_alignment[0]:
        best_alignment = Alignment(
            identity, alignment.get_reference_sequence().upper(), filename)

if best_alignment is not None:
    print(">%s %f" % (best_alignment.filename, best_alignment.indentity))
    print(best_alignment.ref)
else:
    print("NA")
