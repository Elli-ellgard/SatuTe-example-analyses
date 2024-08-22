from Bio import Phylo
from io import StringIO
import pandas as pd
from ete3 import Tree
from collections import Counter
from typing import Tuple, Set
from pandas import DataFrame
import re
import os

def convert_tree_to_newick_string(tree_obj):
    """
    Converts a Tree object back into a Newick string.

    Args:
        tree_obj (Tree): An ETE3 Tree object.

    Returns:
        str: The Newick string representing the tree.
    """
    return tree_obj.write(format=1)

def rename_internal_nodes_pre_order(tree: Tree) -> Tree:
    """
    Modifies the input tree in-place by naming its internal nodes using a pre-order traversal.
    Internal nodes are named as "NodeX*" where X is an incremental number. This function skips renaming
    for nodes that already start with "Node" to avoid redundancy and does not rename numeric node names,
    considering them as not annotated.
    Args:
        tree (Tree): The input tree to be modified.

    Returns:
        Tree: The same tree instance with updated internal node names. Note that this function modifies the tree in-place.
    """
    number_annotated_internal_nodes, number_unannotated_internal_nodes, seen_names = (
        get_set_and_number_annotated_internal_nodes(tree)
    )
    if number_unannotated_internal_nodes != 0:
        for node in tree.traverse("preorder"):
            if node.is_leaf() or node.name.startswith("Node"):
                continue  # Skip leaf nodes and nodes already properly named

            # Rename internal nodes with a placeholder name or numeric name
            if not node.name or node.name.isdigit():
                new_name = f"Node{number_annotated_internal_nodes + 1}*"
                while new_name in seen_names:
                    number_annotated_internal_nodes += 1
                    new_name = f"Node{number_annotated_internal_nodes + 1}*"

                node.name = new_name
                seen_names.add(new_name)
                number_annotated_internal_nodes += 1

    # Optional: Check for any duplicate node names which could cause issues
    check_and_raise_for_duplicate_nodes(tree)

    return tree


def get_set_and_number_annotated_internal_nodes(
    tree: Tree,
) -> Tuple[int, int, Set[str]]:
    """
    Analyzes the internal nodes of a given tree and categorizes them into annotated and unannotated.

    Annotated nodes are those that have names starting with "Node" and are neither empty nor purely numeric.
    Unannotated nodes are either unnamed, numeric named, or have a generic placeholder name.

    The function performs a pre-order traversal of the tree to ensure each node is visited in pre-order sequence.

    Args:
        tree (Tree): The tree to analyze for node annotation status.

    Returns:
        Tuple[int, int, Set[str]]:
            - First int represents the count of annotated internal nodes.
            - Second int represents the count of unannotated internal nodes.
            - Set[str] contains the unique names of the annotated internal nodes.

    This function does not modify the input tree.
    """
    seen_names = set()
    number_annotated = 0
    number_unannotated = 0
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue  # Skip leaf nodes to focus on internal nodes

        if not node.name or node.name.isdigit() or not node.name.startswith("Node"):
            number_unannotated += 1
        else:
            seen_names.add(node.name)
            number_annotated += 1

    return number_annotated, number_unannotated, seen_names


def check_and_raise_for_duplicate_nodes(tree: Tree) -> None:
    """
    Checks the given tree for any duplicate node names and raises an exception if duplicates are found.
    This ensures the tree structure is uniquely identifiable and consistent for downstream analysis.

    Args:
        tree (Tree): The tree to check for duplicate node names.

    Raises:
        ValueError: If a duplicate node name is found in the tree.
    """
    seen_names = set()
    for node in tree.traverse("preorder"):
        if node.name in seen_names:
            raise ValueError(f"Duplicate node name found: '{node.name}'")
        seen_names.add(node.name)


def check_all_internal_nodes_annotated(tree: Tree) -> bool:
    """
    Checks if all internal nodes in the given tree are annotated, where an internal node is considered annotated
    if it has a non-empty, non-numeric name and does not start with a generic prefix like "Node".

    Args:
        tree (Tree): The tree to check for annotated internal nodes.

    Returns:
        bool: True if all internal nodes are annotated, False otherwise.
    """
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue  # Focus on internal nodes
        # An internal node is unannotated if its name is empty, numeric, or a generic placeholder
        if not node.name or node.name.isdigit() or not node.name.startswith("Node"):
            return False
    return True



# import re

# def get_newick_string_from_satute(file_path):
#     """
#     Extracts the Newick string from a .satute file.

#     Args:
#         file_path (str): Path to the .satute file.

#     Returns:
#         str: The Newick string if found, otherwise an empty string.
#     """
#     newick_string = ""

#     # Regular expression pattern to match the Newick tree string
#     newick_pattern = re.compile(r'tree:\s*(\([^;]+;\))', re.DOTALL)

#     try:
#         with open(file_path, 'r') as file:
#             content = file.read()
#             match = newick_pattern.search(content)
#             if match:
#                 newick_string = match.group(1)

#     except FileNotFoundError:
#         print(f"File not found: {file_path}")
#     except Exception as e:
#         print(f"An error occurred: {e}")

#     print(f"new: {newick_string}")

#     return newick_string



def get_newick_string_from_satute(file_path):
    """
    Extracts the Newick string from a .satute file.

    Args:
        file_path (str): Path to the .satute file.

    Returns:
        str: The Newick string if found, otherwise an empty string.
    """
    newick_string = ""

    try:
        with open(file_path, 'r') as file:
            content = file.read()

            # Find the position of "tree:"
            tree_pos = content.find("tree:")
            if tree_pos != -1:
                # Find the start of the Newick string (first '(' after "tree:")
                start_pos = content.find("(", tree_pos)
                # Find the end of the Newick string (first ';' after start_pos)
                end_pos = content.find(";", start_pos)
                
                if start_pos != -1 and end_pos != -1:
                    newick_string = content[start_pos:end_pos+1]

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

    return newick_string


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
    target_node = get_target_node(row["branch"])
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

