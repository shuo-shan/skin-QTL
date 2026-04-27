import json
import os
import pandas as pd

# Directory containing the qc.json files
directory = '/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/QC'

# List to store all QC records
all_records = []

# Iterate through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith("qc.json"):
        file_path = os.path.join(directory, filename)
        with open(file_path, 'r') as file:
            data = json.load(file)
            for record in data:
                # Add a file identifier to each record to trace which file it came from
                record['file'] = filename
                all_records.append(record)

# Flatten the JSON structure and compile into a DataFrame
qc_data = pd.json_normalize(all_records)

# Save to a tab-delimited file
qc_data.to_csv("merged_qc_data.tsv", sep='\t', index=False)

print("QC data has been compiled and saved to 'merged_qc_data.tsv'.")
