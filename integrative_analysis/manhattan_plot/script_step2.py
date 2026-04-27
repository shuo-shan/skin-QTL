import sys
import pandas as pd

# read in file 
# for each unique snp, find the row with the lowest pvalue, and append to output file

# set-up
input_file = '/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/manhattan_plot/temp_result.txt'
output_file = '/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/manhattan_plot/best_associated_PBSeQTL_pairs_and_pval.txt'
# read in file
df = pd.read_csv(input_file, delimiter='\t', header=None, index_col=None)

# 2. for each unique field 1, find the row with the lowest value at field 15, and append to output file
# Group the modeling result by unique SNP values (1st column)
grouped = df.groupby(df.iloc[:,0])

# for each group, find the row with the minimal pvalue (15th column)
result = grouped.apply(lambda group: group[group.iloc[:, 14] == group.iloc[:, 14].min()])
result_focused = result.iloc[:,[0,1,14]]
# write to file
result_focused.to_csv(output_file, sep='\t', header=False, index=False)

# 3. annotate each SNP with its genome location

