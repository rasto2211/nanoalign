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
