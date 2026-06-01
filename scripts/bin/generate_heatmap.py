#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <featurecounts.tsv> <mapping.csv> <gene_list.txt>")
    sys.exit(1)

featurecounts_file = sys.argv[1]
mapping_file = sys.argv[2]
gene_list_file = sys.argv[3]

df = pd.read_csv(featurecounts_file, sep="\t", comment="#")

map_df = pd.read_csv(mapping_file)
map_df = map_df[["ensembl_gene_id", "hgnc_symbol"]].drop_duplicates()

df = df.merge(map_df, left_on="Geneid", right_on="ensembl_gene_id", how="left")
df["GeneName"] = df["hgnc_symbol"].fillna(df["Geneid"])

#TPM Normalization
gene_length = df["Length"] / 1000

expr_cols = [c for c in df.columns if c not in [
    "Geneid", "Chr", "Start", "End", "Strand",
    "Length", "ensembl_gene_id", "hgnc_symbol", "GeneName"
]]

rpk = df[expr_cols].div(gene_length, axis=0)
scaling_factor = rpk.sum(axis=0) / 1e6
tpm = rpk.div(scaling_factor, axis=1)

tpm_df = pd.concat([df[["GeneName"]], tpm], axis=1)

with open(gene_list_file) as f:
    gene_set = set(line.strip() for line in f if line.strip())

tpm_df = tpm_df[tpm_df["GeneName"].isin(gene_set)]

out_prefix = featurecounts_file.replace(".out", "").replace(".tsv", "")

tpm_out = f"{out_prefix}_filtered_TPM_matrix.tsv"
tpm_df.to_csv(tpm_out, sep="\t", index=False)

print(f"TPM matrix saved: {tpm_out}")

#HEATMAP
heatmap_df = tpm_df.set_index("GeneName")

plt.figure(figsize=(12, 8))
plt.imshow(heatmap_df.values, aspect='auto', cmap='viridis_r')

plt.xticks(
    ticks=np.arange(len(heatmap_df.columns)),
    labels=heatmap_df.columns,
    rotation=90
)

plt.yticks(
    ticks=np.arange(len(heatmap_df.index)),
    labels=heatmap_df.index
)

plt.title("Heatmap of TPM normalized gene counts")
plt.xlabel("Samples")
plt.ylabel("Genes")

plt.colorbar(label="TPM")
plt.tight_layout()

heatmap_file = f"{out_prefix}_TPM_heatmap.png"
plt.savefig(heatmap_file, dpi=300)
plt.close()

print(f"Heatmap saved: {heatmap_file}")
