#!/usr/bin/env python3

import pandas as pd
import sys

if len(sys.argv) < 5:
    print(f"Usage: {sys.argv[0]} output.csv mutect2.csv vardict.csv varRNA.csv")
    sys.exit(1)

outfile, mutect_file, vardict_file, varrna_file = sys.argv[1:5]

mutect = pd.read_csv(mutect_file)
vardict = pd.read_csv(vardict_file)
varrna = pd.read_csv(varrna_file)

mutect["caller"] = "Mutect2"
vardict["caller"] = "VarDict"
varrna["caller"] = "VarRNA"

key_cols = ["Chr", "Start", "End", "Ref", "Alt"]

all_variants = pd.concat([mutect, vardict, varrna], ignore_index=True)

grouped = all_variants.groupby(key_cols)
rows = []

for variant, group in grouped:
    row = dict(zip(key_cols, variant))
    callers = []

    for caller in ["Mutect2", "VarDict", "VarRNA"]:
        if caller in group["caller"].values:
            callers.append(caller)

    row["Caller_count"] = len(callers)
    row["Variant_callers"] = ";".join(callers)

    # ❌ SOMATIC_FLAG removed here
    for col in ["FILTER"]:
        vals = []
        for caller in ["Mutect2", "VarDict", "VarRNA"]:
            sub = group[group["caller"] == caller]
            if not sub.empty and col in sub.columns:
                vals.append(str(sub.iloc[0][col]))
        if len(vals) == 1:
            row[col] = vals[0]
        else:
            row[col] = ";".join(vals)

    for col_base in ["REF_COUNT", "ALT_COUNT", "VAF%"]:
        for caller in ["Mutect2", "VarDict", "VarRNA"]:
            colname = f"{col_base}_{caller}"
            sub = group[group["caller"] == caller]
            if not sub.empty and col_base in sub.columns:
                row[colname] = sub.iloc[0][col_base]
            else:
                row[colname] = -1

    for col in ["REF_F2R1", "REF_F1R2", "ALT_F2R1", "ALT_F1R2"]:
        sub = group[group["caller"] == "Mutect2"]
        if not sub.empty and col in sub.columns:
            row[col] = sub.iloc[0][col]
        else:
            row[col] = "-"

    for col in group.columns:
        if col in key_cols + ["caller", "FILTER",
                              "REF_COUNT", "ALT_COUNT", "VAF%",
                              "REF_F2R1", "REF_F1R2", "ALT_F2R1", "ALT_F1R2"]:
            continue
        if col not in row:
            val = group.iloc[0][col]
            try:
                row[col] = pd.to_numeric(val)
            except (ValueError, TypeError):
                row[col] = val

    rows.append(row)

final_df = pd.DataFrame(rows)

final_df.to_csv(outfile, index=False)

print(f"Saved merged file: {outfile}")
