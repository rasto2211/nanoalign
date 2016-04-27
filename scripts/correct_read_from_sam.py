# Takes SAM file with alignment of the read with reference sequence and returns
# the part of ref. seq. that was aligned to the read.
import sys
import pysam
import sam_utils

from collections import namedtuple

file_path = sys.argv[1]
best_alignment, best_identity = sam_utils.get_best_alignment(file_path)

if best_alignment is not None:
    print(">%s %f" % (file_path, best_identity))
    print(best_alignment.get_reference_sequence())
else:
    print("NA")
