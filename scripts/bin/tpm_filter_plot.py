#!/usr/bin/env python

import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# INPUT FILES
# -----------------------------
featurecounts_file = "MYFU_normals_CGH485.tsv"
mapping_file = "/home/hemat/gulnaz/rnaseq_cnv/casper_plot/hg38_casper_base_table.csv"
bed_file = "MECOM.bed"

# -----------------------------
# 1. LOAD FEATURECOUNTS
# -----------------------------
df = pd.read_csv(featurecounts_file, sep="\t", comment="#")

# -----------------------------
# 2. LOAD MAPPING FILE
# -----------------------------
map_df = pd.read_csv(mapping_file)
map_df = map_df[["ensembl_gene_id", "hgnc_symbol"]].drop_duplicates()

df = df.merge(map_df, left_on="Geneid", right_on="ensembl_gene_id", how="left")
df["GeneName"] = df["hgnc_symbol"].fillna(df["Geneid"])

# -----------------------------
# 3. TPM NORMALIZATION (DO THIS BEFORE FILTERING)
# -----------------------------

# Gene length in kb
gene_length = df["Length"] / 1000

# Expression columns (skip metadata columns)
expr_cols = [c for c in df.columns if c not in [
    "Geneid", "Chr", "Start", "End", "Strand",
    "Length", "ensembl_gene_id", "hgnc_symbol", "GeneName"
]]

# RPK
rpk = df[expr_cols].div(gene_length, axis=0)

# Scaling factor per sample
scaling_factor = rpk.sum(axis=0) / 1e6

# TPM
tpm = rpk.div(scaling_factor, axis=1)

# Combine TPM with gene names
tpm_df = pd.concat([df[["GeneName"]], tpm], axis=1)

# -----------------------------
# 4. EXTRACT TARGET GENES FROM BED
# -----------------------------
bed_genes = set()

with open(bed_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 4:
            gene = parts[3].split("_")[0]
            bed_genes.add(gene)

# Filter AFTER TPM
tpm_df = tpm_df[tpm_df["GeneName"].isin(bed_genes)]

# -----------------------------
# 5. SAVE OUTPUT
# -----------------------------
tpm_df.to_csv("26CGH485_filtered_TPM_matrix.tsv", sep="\t", index=False)

print("Done: filtered_TPM_matrix.tsv generated")

# -----------------------------
# 6. HEATMAP GENERATION
# -----------------------------
heatmap_df = tpm_df.set_index("GeneName")

heatmap_log = np.log2(heatmap_df + 1)
data = heatmap_log.values

plt.figure(figsize=(12, 8))
plt.imshow(data, aspect='auto', cmap='viridis_r')

plt.xticks(ticks=np.arange(len(heatmap_df.columns)),
           labels=heatmap_df.columns,
           rotation=90)

plt.yticks(ticks=np.arange(len(heatmap_df.index)),
           labels=heatmap_df.index)

plt.title("TPM Heatmap (log2 transformed)")
plt.xlabel("Samples")
plt.ylabel("Genes")

plt.colorbar(label="log2(TPM + 1)")
plt.tight_layout()

plt.savefig("TPM_heatmap_26CGH485.png", dpi=300)
plt.close()

print("Heatmap saved: TPM_heatmap.png")

# -----------------------------
# 7. HEATMAP (RAW TPM)
# -----------------------------

heatmap_raw = heatmap_df.values

plt.figure(figsize=(12, 8))

plt.imshow(heatmap_raw, aspect='auto', cmap='viridis_r')

plt.xticks(ticks=np.arange(len(heatmap_df.columns)),
           labels=heatmap_df.columns,
           rotation=90)

plt.yticks(ticks=np.arange(len(heatmap_df.index)),
           labels=heatmap_df.index)

plt.title("Heatmap of TPM normalized gene counts")
plt.xlabel("Samples")
plt.ylabel("Genes")

plt.colorbar(label="TPM")

plt.tight_layout()

plt.savefig("TPM_heatmap_raw_26CGH485.png", dpi=300)
plt.close()

print("Heatmap saved: TPM_heatmap_raw.png")
