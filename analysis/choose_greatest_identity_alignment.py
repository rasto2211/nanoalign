import sys

header = input()

greatest_identity = 0
for row in sys.stdin:
    if row == "NA\n":
        print("NA")
        sys.exit(0)
    cells = row.split('\t')
    edit_distance = float(cells[0])
    MID = float(cells[1])
    identity = 1 - edit_distance / MID
    greatest_identity = max(greatest_identity, identity)

print(greatest_identity * 100)
