#! /usr/bin/env python
# This script will convert the metafuse output from hg38 to hg19
import sys
import os
import pandas as pd
import openpyxl
import re
from pyliftover import LiftOver

chainfile = sys.argv[1]	# hg19 to hg38 chain file
metafuse_input = sys.argv[2]	# input excel file
metafuse_output = sys.argv[3]	# output excel file

lo = LiftOver(chainfile)

default_output = ["chrgl", "NA", "NA", "NA"]
def liftover_funtion_ucsc (chromosome, position):
	chrom = str(chromosome)
	pos = int (position)
	#print (chrom, pos)
	output = lo.convert_coordinate(chrom, pos)
	if not output:
		output.append(default_output)
	return output

xl_file = pd.read_excel(metafuse_input, sheet_name=None)
with pd.ExcelWriter(metafuse_output) as writer:
	for sheet_names in xl_file.keys():
		df = pd.read_excel(metafuse_input, sheet_name=sheet_names)
		#for index, row in df.iterrows():
		#	print (index, row['chr1'], row['breakpoint_1'], row['chr2'], row['breakpoint_2'])
		for row in df.index:
			left_chr = df['chr1'][row]
			left_pos = df['breakpoint_1'][row]
			left_convert = liftover_funtion_ucsc(left_chr, left_pos)

			right_chr = df['chr2'][row]
			right_pos = df['breakpoint_2'][row]
			right_convert = liftover_funtion_ucsc(right_chr, right_pos)

			df.at[row, 'chr1'] = left_convert[0][0]
			df.at[row, 'breakpoint_1']= left_convert[0][1]
			df.at[row, 'chr2'] = right_convert[0][0]
			df.at[row, 'breakpoint_2'] = right_convert[0][1]

		df.rename(columns = {'exon1':'exon1_hg19', 'exon2':'exon2_hg19'}, inplace = True)
		df.to_excel(writer, sheet_name=sheet_names, index=False)	
