# Samples reads from CSV file produced by `merge_csv.py` and
# `csv_with_read_stats.py`. First argument is the number of samples and 2nd is
# the path to CSV file.

import sys
import random

import numpy as np
import pandas as pd

n_samples = int(sys.argv[1])
input_file = sys.argv[2]

csv = pd.read_csv(input_file, na_values=["NA"])
csv = csv.rename(columns=lambda x: x.strip())

csv["match"] = csv["M"] - (csv["edit_dist"] - csv["I"] - csv["D"])
csv["id"] = csv["match"] / (csv["edit_dist"] + csv["match"])
csv["id_clip"] = csv["match"] / \
    (csv["edit_dist"] + csv["match"] + csv["S_CLIP"] + csv["H_CLIP"])

selected_paths = np.random.choice(csv.merged_path.unique(),
                                  n_samples,
                                  replace=False)
selected_rows = csv[csv.merged_path.isin(selected_paths)]

grouped_id_clip = selected_rows[
    ["merged_path", "id_clip"]].groupby('merged_path')
grouped_id = selected_rows[["merged_path", "id"]].groupby('merged_path')

id_clip_df = pd.DataFrame()
for name, group in grouped_id_clip:
        id_clip_df[name] = list(group.id_clip)
id_clip_df.to_csv('id_clip_samples.csv', index=False)

id_df = pd.DataFrame()
for name, group in grouped_id:
        id_df[name] = list(group.id)
id_df.to_csv('id_samples.csv', index=False)
