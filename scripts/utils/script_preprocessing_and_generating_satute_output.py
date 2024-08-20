import os
import subprocess
import shutil
import time

import sys
# Add the directory containing the satute module to the Python path
sys.path.append('/home/elgert/Desktop/Cassius/version_2024_08_19')
# Now you can import satute.cli
import satute.cli


from  utils.script_analyses_utils import (
    find_file_with_suffix_in_directory
)


def run_satute(args):
    """
    Runs the satute CLI with the provided arguments, ensuring the 'quiet' 
    argument is set, and adds a small delay to ensure files are written.

    Args:
        args (list or None): Command-line arguments to pass to satute.cli.main. If None, an empty list is used.

    Returns:
        bool: The result of satute.cli.main(args).
    """
    if args is None:
        args = []

    if not isinstance(args, list):
        raise TypeError("args must be a list or None")

    # Ensure '--quiet' is in the arguments
    if '-quiet' not in args:
        args.append('-quiet')

    # Add a small delay to ensure files are written
    time.sleep(2)

    # Call the main function of satute.cli with the arguments
    return satute.cli.main(args)


def  run_satute_for_directory(folder_path, path_iqtree, path_python, path_satute, alpha):

    print(folder_path)

    run_satute(
        [
            "-dir",
            folder_path,
            "-alpha",
            alpha,
        ]
    )

    # arguments = [
    #     path_python, 
    #     path_satute,
    #     "-dir",
    #     folder_path,
    #     "-alpha",
    #     alpha,
    # ]
    # print(arguments)
    # run_external_command(arguments)
        

        



# def preprocessing_and_generating_satute_output(data_name, fasta_file, output_dir, path_iqtree, path_python, path_satute, alpha):
#     # Print some information
#     print("Considered file: ", data_name)
#     print("")

#     # Generate the directory for analysis
#     current_dir= generate_directory("full_genome", fasta_file, output_dir)

#     # Run Satute
#     run_satute("genome", current_dir, path_iqtree, path_python, path_satute, alpha)

#     return current_directory




""" Preprocessing for branch specific analyses """

def generate_edge_directories(edge_list, source_folder, output_dir):
    # Iterate over each gene name in the list
    for edge_name in edge_list:
        print(f"Generate directory for {edge_name.upper()}")

        # Create a directory for the gene
        edge_directory = os.path.join(output_dir, edge_name)
        os.makedirs(edge_directory, exist_ok=True)

        # Find the corresponding FASTA file in the source folder
        fasta_file = find_file_with_suffix_in_directory(".fasta", source_folder)
        
        if fasta_file:
            # Copy the FASTA file to the gene directory
            shutil.copy(fasta_file, edge_directory)

        else: 
            print(f"Fasta file does not exists!")

def  run_satute_for_edge(edge_name, tree_file, folder_path, path_iqtree, alpha,model):
    # Run Satute using  alignment and tree
    log_file = find_file_with_suffix_in_directory("satute.log", folder_path)
    if not log_file:
        fasta_file_aln = find_file_with_suffix_in_directory(".fasta", folder_path)
        arguments = [
                "-iqtree",
                path_iqtree,
                "-msa",
                fasta_file_aln,
                "-tree",
                tree_file,
                "-model",
                model, 
                "-alpha",
                alpha,
                "-edge",
                f"{edge_name}"
        ]
        run_satute(arguments)
        
    else:
        print("SatuTe run is already there!")
        


def branch_specific_preprocessing(input_dir, edges_dict, output_dir, path_iqtree, alpha, model):
    # Print some information
    print("Number of considered egdes: ", len(edges_dict))
    print("Considered edges: ", edges_dict.keys())
    print("")


    # Generate the directory and complementary alignments for all genes 
    generate_edge_directories(edges_dict.keys(), input_dir, output_dir)
    print("")

    tree_file = find_file_with_suffix_in_directory(".treefile",input_dir) 

    # Generate necessary data 
    for edge_folder in edges_dict.keys():
        print(f"Run SatuTe for {edge_folder.upper()}")
        folder_path = os.path.join(output_dir, edge_folder)
        run_satute_for_edge(edges_dict[edge_folder], tree_file, folder_path, path_iqtree, alpha, model)
        print("")