import requests, json,re, os, sys,subprocess

##################################
def fetch(url):
  headers = {'accept': 'application/json'}
  response = requests.get(url, headers=headers)
  data = response.json()
  return(data)

## get the URL page
#url="https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=IMR-90"
url="https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=TF+ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.cell_slims=epithelial+cell&biosample_ontology.cell_slims=fibroblast&biosample_ontology.cell_slims=connective+tissue+cell&biosample_ontology.cell_slims=keratinocyte"
data=fetch(url+"&limit=all&format=json")

f=open('metadata.txt','w')

for exp in data.get('@graph', []):  # Use .get() with a default empty list for safety
    exp_ID = exp['accession']  # Use .get() with a default empty string
    exp_url = f"https://www.encodeproject.org/experiments/{exp_ID}/?format=json"
    headers = {'accept': 'application/json'}
    exp_response = requests.get(exp_url, headers=headers)
    exp_metadata = exp_response.json()
    for entry in exp_metadata.get("files", []):  # Use .get() with a default empty list
        # Get peak bed file
        if entry.get("file_type") == "bed narrowPeak":
            # Initialize a flag to False, it will help to check conditions progressively
            encode4_condition_met = False
            # Check 'analyses' key existence and non-emptiness, then further check within
            if 'analyses' in entry and entry['analyses'] and 'pipeline_award_rfas' in entry['analyses'][0]:
                if entry['analyses'][0]["pipeline_award_rfas"]:  # Checks list is not empty
                    encode4_condition_met = entry['analyses'][0]["pipeline_award_rfas"][0] == 'ENCODE4'
            # Proceed if ENCODE4 condition met
            if encode4_condition_met:
                # Get the preferred default file
                if entry.get("preferred_default") == True:
                    if entry.get("status") == "released":
                        url = f"https://www.encodeproject.org/files/{entry['accession']}/@@download/{entry['accession']}.bed.gz\n"
                        target = entry.get('target', 'No Target Specified')  # Example of using .get() with a different default
                        biosample_summary = exp.get("biosample_summary", "No Biosample Summary")  # Default if not present
                        assembly = entry.get("assembly", "No Assembly")  # Default if not present
                        simple_biosample_summary = entry.get('simple_biosample_summary', 'No Simple Biosample Summary')  # Default if not present
                        f.write(f"{exp_ID};{entry['accession']};{target};{biosample_summary};{assembly};{url};{simple_biosample_summary}")


#for exp in data['@graph']:
#  exp_ID=exp['accession']
#  exp_url = "https://www.encodeproject.org/experiments/" + exp_ID + "/?format=json"
#  headers = {'accept':'application/json'}
#  exp_response = requests.get(exp_url, headers = headers)
#  exp_metadata = exp_response.json()
#  for entry in exp_metadata["files"]:
#    #get peak bed file
#    if (entry["file_type"] == "bed narrowPeak") :
#      # get encode 4 annotation
#      if 'analyses' in entry and "pipeline_award_rfas" in entry['analyses'][0]:
#          if len(entry['analyses'][0]["pipeline_award_rfas"])>0 and entry['analyses'][0]["pipeline_award_rfas"][0] == 'ENCODE4':
#              # get the preffered default  file
#              if "preferred_default" in entry and entry["preferred_default"] == True:
#                  if entry["status"] == "released":
#                      url = "https://www.encodeproject.org/files/" + entry['accession'] + "/@@download/" + entry['accession'] + ".bed.gz" + '\n'
#                      f.write(exp['accession'] + ';' + entry['accession'] + ';' + entry['target'] + ';' + exp["biosample_summary"] + ';' + entry["assembly"] + ';' + url + ';' + entry['simple_biosample_summary'])

f.close()
