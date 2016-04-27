import sys
import random

def remove_percent(token):
    if token.endswith("%"):
        return token[:-1]

    return token

n_samples = int(sys.argv[1])

sample_paths = random.sample([path.strip() for path in sys.stdin], n_samples)

cols = []
for path in sample_paths:
    file = open(path)
    header = file.readline()
    indentities = [line.strip() for line in file]
    cols.append(indentities)

# Transpose matrix of cols.
csv_matrix = list(map(list, zip(*cols)))

# Print header.
print(",".join(["r" + str(i + 1) for i in range(n_samples)]))

for row in csv_matrix:
    print(",".join([remove_percent(cell) for cell in row]))
