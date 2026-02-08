import os, glob, sys
import pandas as pd
from pathlib import Path

def find_gather_csv_files(base_dir, check_table):
    samples = set()
    if check_table != "None" and check_table:
        with open(check_table, "r") as f:
            for line in f:
                if line.startswith("#"): continue
                tokens = line.strip().split()
                samples.add(tokens[1])                

    base_path = Path(base_dir)
    csv_files = []

    for subdir in base_path.iterdir():
        if subdir.is_dir() and subdir.name.startswith("gather"):
            for csv_file in subdir.rglob("*.csv"):
                accession = Path(csv_file).stem
                if accession not in samples: continue
                csv_files.append(csv_file)
    return csv_files

# Receive inputs from Snakemake
output_path = sys.argv[2]
sample_paths = find_gather_csv_files(sys.argv[1], sys.argv[3])

def get_relab(sample_path):
    sample_file = os.path.basename(sample_path)
    sample = os.path.splitext(sample_file)[0]
    df = pd.read_csv(sample_path)
    if "potential_false_negative" in df.columns:
        df = df[df.potential_false_negative == False]
    df = df[df.intersect_bp >= 25000]
    if df.empty:
        return pd.DataFrame()
    df[sample] = df["f_unique_weighted"] / df["f_unique_weighted"].sum()
    df = df[["name", sample]].set_index("name")
    df = df.T.groupby(axis=1, level=0).sum().T
    return df


# Aggregate all dataframes
all_dfs = []
for path in sample_paths:
    try:
        df = get_relab(path)
        if not df.empty:
            all_dfs.append(df)
    except Exception as e:
        print(f"Failed to process {path}: {e}", file=sys.stderr)

if all_dfs:
    final_df = pd.concat(all_dfs, axis=1).fillna(0)
    final_df.index = final_df.index.astype(str).str.rjust(10, "0")
    final_df.to_csv(output_path)
else:
    print("No valid dataframes to concatenate.", file=sys.stderr)
    pd.DataFrame().to_csv(output_path)
