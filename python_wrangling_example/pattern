#!/usr/bin/env python3

# this is a pattern that one could follow to compare vcfs and a standard

import numpy as np
import pandas as pd

minimap2 = pd.read_csv('processed.minimap2_nanopore.vcf.csv', sep=",", index_col="PROV")
ngmlr = pd.read_csv('processed.ngmlr_nanopore.vcf.csv', sep=",", index_col="PROV")
winnowmap2 = pd.read_csv('processed.winnowmap_nanopore.vcf.csv', sep=",", index_col="PROV")

merged_np = pd.concat([minimap2,ngmlr,winnowmap2])

bins = [50, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000]
merged_np['BINS'] = pd.cut(merged_np['ABS_SVLEN'], bins)

np_del = merged_np[merged_np["SVTYPE"] == "DEL"]
np_ins = merged_np[merged_np["SVTYPE"] == "INS"]
np_inv = merged_np[merged_np["SVTYPE"] == "INV"]
np_dup = merged_np[merged_np["SVTYPE"] == "DUP"]

new_truth = pd.read_csv('new_truth.csv',sep=',', index_col="PROV")

del_truth = new_truth[new_truth["SVTYPE"] == "deletion"]
ins_truth = new_truth[new_truth["SVTYPE"] == "insertion"]
inv_truth = new_truth[new_truth["SVTYPE"] == "inversion"]
dup_truth = new_truth[new_truth["SVTYPE"] == "duplication"]

all_del = pd.concat([del_truth,np_del])
all_del.to_csv('np_truth_del_binned.csv')

all_ins = pd.concat([ins_truth,np_ins])
all_ins.to_csv('np_truth_ins_binned.csv')

all_inv = pd.concat([inv_truth,np_inv])
all_inv.to_csv('np_truth_inv_binned.csv')

all_dup = pd.concat([dup_truth,np_dup])
all_dup.to_csv('np_truth_dup_binned.csv')
