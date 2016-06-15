# Merges multiple CSV files of the same type. The first column of the resulting
# CSV contains name of the file from which the row comes from. All the other
# columns are copied from the original CSV files. The list of files to merge is
# read from stdin.

import sys

first_file = True
for line in sys.stdin:
    path = line.strip()
    file = open(path)
    header = "merged_path, " + file.readline()

    if first_file:
        print(header.strip())
        first_file = False

    for file_line in file:
        print("%s,%s" %(path, file_line.strip()))
