from Bio import Phylo
from io import StringIO
import pandas as pd
from pandas import DataFrame
import re
import os

# we  use the known structure of the nexus file ... everything else was not working 
def extract_newick_from_nexus(nexus_file):
    with open(nexus_file, 'r') as nexus_handle:
        # Read the Nexus file content
        nexus_content = nexus_handle.read()

        # Find the start and end indices of the tree section
        start_index = nexus_content.find("BEGIN TREES;")
        end_index = nexus_content.find("END;", start_index)

        # Extract the tree section
        tree_section = nexus_content[start_index:end_index]

        # Extract the Newick string from the tree section
        newick_start = tree_section.find("(")
        newick_end = tree_section.find(";", newick_start) + 2  # Include the closing bracket and semicolon
        newick_string = tree_section[newick_start:newick_end]
        
        return newick_string


            
def get_plain_tree(directory_path):
    # Get a list of nexus files
    file_suffix = ".satute.nex"
    list_nexus_files = [f for f in os.listdir(directory_path) if f.endswith(file_suffix)]
    nexus_file = os.path.join(directory_path, list_nexus_files[0])
    print(nexus_file)
    newick_string = extract_newick_from_nexus(nexus_file)
    cleaned_newick = remove_metadata_from_newick(newick_string)

    # print("Original Newick string:")
    print(newick_string)

    # print("\nNewick string without metadata:")
    # print(cleaned_newick)
    return cleaned_newick

def list_files_with_suffix(folder, suffix):
    matching_files = [f for f in os.listdir(folder) if f.endswith(suffix)]
    return matching_files


def get_tree_file(directory_path):
    # Get a list of files
    file_suffix = "satute_tree.tree"
    list_nexus_files = list_files_with_suffix(directory_path, file_suffix)
    if not list_nexus_files:
        file_suffix = "satute.nex"
        list_nexus_files = list_files_with_suffix(directory_path, file_suffix)
    return os.path.join(directory_path, list_nexus_files[0])
    

""" Extract plain newick string of considered tree """
def is_nexus(content):
    return bool(re.search(r'#NEXUS', content, re.IGNORECASE))

def extract_newick(content):
    newick_match = re.search(r'\s*([^#].*;)\s*', content, re.DOTALL)
    if newick_match:
        return newick_match.group(1)
    return None

def remove_metadata_from_newick(newick_string):
    # Remove metadata (information of the form [& ...])
    cleaned_newick = ''
    inside_metadata = False

    for char in newick_string:
        if char == '[':
            inside_metadata = True
        elif char == ']':
            inside_metadata = False
        elif not inside_metadata:
            cleaned_newick += char

    return cleaned_newick

def get_newick_string(file_path):
    with open(file_path, 'r') as file:
        file_content = file.read()

    # Test if it's a Nexus file
    if is_nexus(file_content):
        #print(f'The file "{file_path}" is a Nexus file.')
        # Extract Newick string
        newick_string = extract_newick_from_nexus(file_path)
    else:
        # Assume the whole content is a Newick string
        newick_string = file_content
        
    if newick_string:
        cleaned_newick = remove_metadata_from_newick(newick_string)
        return(cleaned_newick)
    else:
        # No newick string in the file
        return("")
    

""" Write new metadata to newick string"""

def get_target_node(edge):
    return edge.split(",")[0].replace("(", "").strip()


def map_values_to_newick(newick: str, results_per_branch: DataFrame) -> str:
    """
    Maps values from the DataFrame onto the Newick string by dynamically modifying the metadata for each node.

    Args:
    - newick (str): The Newick string to be modified.
    - results_per_branch (DataFrame): DataFrame containing the data to map onto the Newick string.

    Returns:
    - str: Modified Newick string with updated metadata.
    """
    for index, row in results_per_branch.iterrows():
        newick = update_node_metadata(newick, row, results_per_branch.columns)

    return newick


def update_node_metadata(newick: str, row: DataFrame, columns: list) -> str:
    """
    Updates the metadata for a single node in the Newick string based on the DataFrame row.

    Args:
    - newick (str): Current Newick string.
    - row (DataFrame): A single row from the DataFrame.
    - columns (list): List of columns in the DataFrame.

    Returns:
    - str: Newick string with updated metadata for the node.
    """
    target_node = get_target_node(row["edge"])
    escaped_target_node = re.escape(target_node)
    meta_data = create_meta_data_string(row, columns)
   
    newick = insert_metadata_into_newick(newick, escaped_target_node, meta_data)
    return newick


def create_meta_data_string(row: DataFrame, columns: list) -> str:
    """
    Creates a metadata string from a DataFrame row.

    Args:
    - row (DataFrame): A single row from the DataFrame.
    - columns (list): List of columns in the DataFrame.

    Returns:
    - str: Metadata string.
    """
    meta_data_parts = [f"{col}={row[col]}" for col in columns if col != "edge"]
    return "&" + ",".join(meta_data_parts) if meta_data_parts else ""


def insert_metadata_into_newick(newick: str, target_node: str, meta_data: str) -> str:
    """
    Inserts metadata into the Newick string for a specific node.

    Args:
    - newick (str): The Newick string.
    - target_node (str): Target node for which metadata needs to be updated.
    - meta_data (str): Metadata string to be inserted.

    Returns:
    - str: Newick string with updated metadata.
    """
    pattern_with_brackets = re.compile(
        rf"({target_node}:\d+(\.\d+)?(e-?\d+)?)\[([^\]]+)\]"
    )
    pattern_without_brackets = re.compile(rf"({target_node}:\d+(\.\d+)?(e-?\d+)?)")

    if pattern_with_brackets.search(newick):
        return pattern_with_brackets.sub(rf"\1[\4,{meta_data}]", newick)
    else:
        return pattern_without_brackets.sub(rf"\1[{meta_data}]", newick)


def write_nexus_file(
    newick_string: str, file_name: str, results_data_frame: DataFrame
):
    """
    Writes the given Newick string as a Nexus file after mapping values.
    Args:
    - newick_string (str): Newick formatted string representing the tree.
    - file_name (str): Name of the file to write.
    - results_data_frame (DataFrame): DataFrame with the data to map onto the Newick string.
    """
    try:
        # Map values to Newick string
        mapped_newick_string = map_values_to_newick(newick_string, results_data_frame)

        # Write the Nexus file
        with open(file_name, "w") as nexus_file:
            nexus_file.write("#NEXUS\n")
            nexus_file.write("BEGIN TREES;\n")
            nexus_file.write(f"Tree tree1 = {mapped_newick_string}\n")
            nexus_file.write("END TREES;\n")
    except Exception as e:
        print(f"Error writing Nexus file: {e}")

