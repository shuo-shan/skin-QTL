import json
import os
import pandas as pd

# Directory containing the ancestry.json files
directory = '/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/072024_donors/ancestry'

# List to store all ancestry records
all_records = []

# Iterate through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith("ancestry.json"):
        file_path = os.path.join(directory, filename)
        with open(file_path, 'r') as file:
            data = json.load(file)
            record = {
                'file': filename,
                'ancestry_metadata_id': data.get('ancestry_metadata_id', None)
            }
            # Add all ancestry and ancestry_raw values
            record.update(data.get('ancestry', {}))
            record.update(data.get('ancestry_raw', {}))
            all_records.append(record)

# Compile all records into a DataFrame
ancestry_data = pd.DataFrame(all_records)

# Save to a tab-delimited file, filling missing values with 'NA'
ancestry_data.to_csv("merged_ancestry_data.tsv", sep='\t', index=False, na_rep='NA')

print("Ancestry data has been compiled and saved to 'merged_ancestry_data.tsv'.")

