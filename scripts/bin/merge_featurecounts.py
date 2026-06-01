#!/usr/bin/env python3

import sys
import os
import re

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <PON_featurecounts> <CASE_featurecounts> <output>")
    sys.exit(1)

pon_file = sys.argv[1]
case_file = sys.argv[2]
out_file = sys.argv[3]

def clean_name(raw_name):
    name = os.path.basename(raw_name)
    name = re.sub(r'\.sorted\.bam$|\.bam$', '', name)
    return name

#Read PON
pon_data = {}
pon_header = None

with open(pon_file) as f:
    lines = f.readlines()

pon_header = lines[1].strip().split("\t")

for line in lines[2:]:
    parts = line.strip().split("\t")
    gene = parts[0]
    pon_data[gene] = parts

#Read CASE
case_data = {}
case_name = None

with open(case_file) as f:
    lines = f.readlines()

case_header = lines[1].strip().split("\t")

raw_case_name = case_header[6]
case_name = clean_name(raw_case_name)

for line in lines[2:]:
    parts = line.strip().split("\t")
    gene = parts[0]
    count = parts[6]
    case_data[gene] = count

with open(out_file, "w") as out:

    new_header = pon_header[:6] + [case_name] + pon_header[6:]
    out.write("\t".join(new_header) + "\n")

    missing = 0
    for gene, pon_row in pon_data.items():
        if gene in case_data:
            case_count = case_data[gene]
            new_row = pon_row[:6] + [case_count] + pon_row[6:]
            out.write("\t".join(new_row) + "\n")
        else:
            missing += 1

if missing > 0:
    print(f"Warning: {missing} genes present in PON but missing in case file", file=sys.stderr)

print(f"Merged file written to: {out_file}")
