# Merges multiple CSV files of the same type. The first column of the resulting
# CSV contains name of the read from which the row comes from. 
# All the other # columns are copied from the original CSV files. 
# The list of files to merge is read from stdin.

# TODO merge_csv.py and this file very similar. Merge them into one script.

import sys
import os

first_file = True
for line in sys.stdin:
    file_path = line.strip()
    read_name = os.path.basename(file_path)[:-len("_intersection.csv")]
    opened_file = open(file_path)

    # Print header only once.
    header = opened_file.readline().strip()
    if first_file:
        print("read_name,%s" %header)
        first_file = False

    for csv_line in opened_file:
        print("%s,%s" % (read_name, csv_line.strip()))
