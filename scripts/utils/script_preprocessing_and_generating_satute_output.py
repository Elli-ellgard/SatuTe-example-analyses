import os
import shutil

from  utils.script_analyses_utils import (
    find_file_with_suffix_in_directory,
    run_external_command
)


""" Preprocessing for branch specific analyses """

def generate_edge_directories(edge_list, output_dir, source_folder=None, fasta_file=None):
    # Iterate over each gene name in the list
    for edge_name in edge_list:
        print(f"Generate directory for {edge_name.upper()}")

        # Create a directory for the gene
        edge_directory = os.path.join(output_dir, edge_name)
        os.makedirs(edge_directory, exist_ok=True)

        # Determine the FASTA file to copy
        if fasta_file:
            # Directly use the provided FASTA file
            file_to_copy = fasta_file
        elif source_folder:
            # Find the corresponding FASTA file in the source folder
            file_to_copy = find_file_with_suffix_in_directory(".fasta", source_folder)
        else:
            print(f"No source folder or FASTA file provided for {edge_name}!")
            continue

        if file_to_copy and os.path.isfile(file_to_copy):
            # Copy the FASTA file to the gene directory
            shutil.copy(file_to_copy, edge_directory)
        else:
            print(f"FASTA file does not exist for {edge_name}!")

def  run_satute_for_edge(edge_name, tree_file, folder_path, path_iqtree, alpha,model):
    # Run Satute using  alignment and tree
    log_file = find_file_with_suffix_in_directory("satute.log", folder_path)
    if not log_file:
        fasta_file_aln = find_file_with_suffix_in_directory(".fasta", folder_path)
        arguments = [
                "satute",
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
        run_external_command(arguments)
        
    else:
        print("SatuTe run is already there!")
     
def branch_specific_preprocessing(input_dir=None, fasta_file=None, tree_file=None, edges_dict=None, output_dir=None, path_iqtree=None, alpha=None, model=None):
    # Validate input: either input_dir or both fasta_file and tree_file must be provided
    if not input_dir and not (fasta_file and tree_file):
        raise ValueError("You must provide either input_dir or both fasta_file and tree_file.")
    
    if not edges_dict or not output_dir or not path_iqtree or alpha is None or model is None:
        raise ValueError("Required parameters are missing.")
    
    # Print some information
    print("Number of considered edges: ", len(edges_dict))
    print("Considered edges: ", edges_dict.keys())
    print("")

    # Determine tree file if not provided but input_dir is
    if not tree_file and input_dir:
        tree_file = find_file_with_suffix_in_directory(".treefile", input_dir)
        if not tree_file:
            raise ValueError("Tree file not found in the input directory.")
    
    # Generate the directory and complementary alignments for all genes 
    generate_edge_directories(edges_dict.keys(), output_dir, source_folder=input_dir, fasta_file=fasta_file)
    print("")

    # Generate necessary data 
    for edge_folder in edges_dict.keys():
        print(f"Run SatuTe for {edge_folder.upper()}")
        folder_path = os.path.join(output_dir, edge_folder)
        run_satute_for_edge(edges_dict[edge_folder], tree_file, folder_path, path_iqtree, alpha, model)
        print("")