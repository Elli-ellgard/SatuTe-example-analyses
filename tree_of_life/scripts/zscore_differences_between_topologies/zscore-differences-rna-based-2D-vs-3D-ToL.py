import os
import sys
# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../scripts/')))


from utils.script_analyses_utils import (
    find_file_with_suffix, 
    find_files_with_suffix_in_directory
)

from utils.script_preprocessing_and_generating_satute_output import (
    branch_specific_preprocessing,
)

from zscore_differences_between_topologies.script_per_region_zscore_differences_between_topologies import (
    per_region_analysis,
    per_region_zscores_toplogies,
)


if __name__ == "__main__":

    # Set paths to IQ-TREE and Python executable
    path_iqtree = "iqtree2"

    # Set significance level
    alpha = str(0.01)
    
    # Get the current working directory
    current_directory = os.path.join(os.getcwd(), "../../")

    """" Preprocessing and Generation of SatuTe output files"""

    # Used evolutionray model
    model = "GTR+G4"

    folder_dict = {}

    # original tree analysis / input files
    dataname1 = "rRNA_based_3D_tree"
    input_dir = os.path.join(current_directory, "data", dataname1)
    region_annotation_file = os.path.join(input_dir, "rRNA_region_annotation.csv")
    edges_dict1= {
        "branch_to_eukaryota": "(Node1612*, Node77*)",
    }

    # Preprocessing: generate directories + run SatuTe and IQ-Tree
    satute_dir_topology1= os.path.join(current_directory,"SatuTe_results", dataname1)
    os.makedirs(satute_dir_topology1, exist_ok=True) 

    branch_specific_preprocessing(
        input_dir=input_dir, 
        edges_dict=edges_dict1, 
        output_dir=satute_dir_topology1, 
        path_iqtree=path_iqtree, 
        alpha=alpha, 
        model =model)

    folder_dict[dataname1] =  os.path.join(satute_dir_topology1, "branch_to_eukaryota")
   
    # rearranged  tree
    dataname2 = "rearranged_rRNA_based_2D_tree"
    fasta_file = find_file_with_suffix("", ".fasta", input_dir)
    rearranged_tree_file = os.path.join(current_directory, "data/rearrangedTrees/16SRNAhugsEtAl-SSU-MSA-Suppl3-2D-joined_with_Loki.tree")
    edges_dict2= {
        "branch_to_eukaryota": "(NodeInserted, Node1710*)",
    }
    satute_dir_topology2= os.path.join(current_directory,"SatuTe_results", dataname2)
    os.makedirs(satute_dir_topology2, exist_ok=True) 

    branch_specific_preprocessing(
        fasta_file=fasta_file,
        tree_file=rearranged_tree_file,
        edges_dict=edges_dict2, 
        output_dir=satute_dir_topology2, 
        path_iqtree=path_iqtree, 
        alpha=alpha, 
        model =model)
    
    folder_dict[dataname2] =  os.path.join(satute_dir_topology2, "branch_to_eukaryota")

    ''' Comparison region zscores '''

    # Specify the path to output
    output_dir = os.path.join(current_directory,"results_zscore_differences_topologies/", dataname1)
    os.makedirs(output_dir, exist_ok=True)
 
    # Results dictionary
    region_zscores = {}

    for dataname, folder in folder_dict.items():
        region_zscores[dataname] = per_region_analysis(
            annotation_file=region_annotation_file, 
            satute_input_dir=folder, 
            results_dir=output_dir, 
            data_name=dataname
        )

    per_region_zscores_toplogies(region_zscores, output_dir, title="rna_based_3D_minus_2D_Eukaryota_branch")

