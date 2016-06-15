# Takes SAM file with alignment of the read to the reference sequence and returns
# the part of ref. seq. that was aligned to the read.

import sys
import pysam
import sam_utils

from collections import namedtuple

file_path = sys.argv[1]
best_alignment = sam_utils.get_best_alignment(file_path)

print(">%s" %file_path)
if best_alignment is not None:
    print(best_alignment.get_reference_sequence().upper())
else:
    print("NA")
