#python3
import dropbox

def list_folder_recursive(dbx, folder_path, output_file, shared_links_map):
    try:
        # Initial call to list files in the folder
        response = dbx.files_list_folder(folder_path)

        # Process entries
        process_entries(dbx, response.entries, output_file, shared_links_map)

        # If there are more entries, continue listing and processing
        while response.has_more:
            response = dbx.files_list_folder_continue(response.cursor)
            process_entries(dbx, response.entries, output_file, shared_links_map)
    except Exception as e:
        print("Error:", e)

def process_entries(dbx, entries, output_file, shared_links_map):
    for entry in entries:
        # Check if the entry is a file
        if isinstance(entry, dropbox.files.FileMetadata):
            shared_link_url = get_or_create_shared_link(dbx, entry, shared_links_map)
            # Write the line to the file
            file_name = entry.path_display.split('/')[-1]
            output_file.write(f"{file_name}\t{entry.path_display}\t{shared_link_url}\n")
        elif isinstance(entry, dropbox.files.FolderMetadata):
            # Recursively list this folder
            list_folder_recursive(dbx, entry.path_lower, output_file, shared_links_map)

def get_or_create_shared_link(dbx, entry, shared_links_map):
    # Check the cache first
    if entry.path_lower in shared_links_map:
        return shared_links_map[entry.path_lower]
    else:
        # Check if a shared link already exists for the file
        links = dbx.sharing_list_shared_links(path=entry.path_lower).links
        if len(links) == 0:
            # No shared link exists; create a new one
            shared_link_metadata = dbx.sharing_create_shared_link_with_settings(entry.path_lower)
            shared_link_url = shared_link_metadata.url
        else:
            # Use the existing shared link
            shared_link_url = links[0].url  # Assuming using the first link if multiple
        
        # Update the cache
        shared_links_map[entry.path_lower] = shared_link_url
        return shared_link_url

# Initialize a Dropbox object using your access token
token = os.environ["DROPBOX_TOKEN"] 
dbx = dropbox.Dropbox(token)

# Specify your Dropbox folder path here, use empty string for root
folder_path = ''

# Open a file for writing
with open('plot_url_lookup.txt', 'w') as output_file:
    # Dictionary to cache shared link URLs and avoid repeated API calls for the same item
    shared_links_map = {}
    # Start recursive listing
    list_folder_recursive(dbx, folder_path, output_file, shared_links_map)











#
#import dropbox
#
## Initialize a Dropbox object using your access token
#dbx = dropbox.Dropbox(token)
#
## Specify your Dropbox folder path here
#folder_path = ''
#
#with open('plot_url_lookup.txt', 'w') as output_file:
#    try:
#        # List files in the folder
#        for entry in dbx.files_list_folder(folder_path).entries:
#            # Check if the entry is a file
#            if isinstance(entry, dropbox.files.FileMetadata):
#                # Check if a shared link already exists for the file
#                links = dbx.sharing_list_shared_links(path=entry.path_lower).links
#                shared_link_url = ""
#                if len(links) == 0:
#                    # No shared link exists; create a new one
#                    shared_link_metadata = dbx.sharing_create_shared_link_with_settings(entry.path_lower)
#                    shared_link_url = shared_link_metadata.url
#                else:
#                    # Use the existing shared link
#                    shared_link_url = links[0].url  # Assuming using the first link if multiple
#
#                # Write the line to the file
#                output_file.write(f"{entry.name}\t{shared_link_url}\n")
#    except Exception as e:
#        print("Error:", e)
