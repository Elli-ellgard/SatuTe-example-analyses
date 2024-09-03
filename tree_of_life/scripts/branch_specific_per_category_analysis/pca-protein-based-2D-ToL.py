import os
import sys
# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../scripts/')))


from utils.script_analyses_utils import (
    find_file_with_suffix, 
)

from utils.script_preprocessing_and_generating_satute_output import (
    branch_specific_preprocessing,
)

from per_category_analysis.script_category_analysis import (
    summary_category_zscore_per_branch,
)

  

if __name__ == "__main__":

    # Set paths to IQ-TREE and Python executable
    path_iqtree = "iqtree2"

    # Set significance level
    alpha = str(0.01)

    """ Set directories and input files """
    # Get the current working directory
    current_directory = os.path.join(os.getcwd(), "../../")


    # dictionary of considered edges
    edges_dict= {
        "branch_to_eukaryota": "(Node2691*, Node2692*)",
        "branch_to_yeast": "(Eukaryota_Opisthokonta_Fungi_Dikarya_Ascomycota_saccharomyceta_Saccharomycotina_Saccharomycetes_Saccharomycetales_Saccharomycetaceae_Saccharomyces_cerevisiae, Node2721*)",
    }

    # Specify the path to your considered data 
    data_name = "protein_based_2D_tree"
    input_dir = os.path.join(current_directory, "data", data_name)
    fasta_file = find_file_with_suffix(data_name, ".fasta", input_dir)
    model = "LG+G4"
    

    satute_output_dir = os.path.join(current_directory,"SatuTe_results", data_name)
    os.makedirs(satute_output_dir, exist_ok=True)

    # Preprocessing: generate directories + run SatuTe and IQ-Tree
    branch_specific_preprocessing(
        input_dir=input_dir, 
        edges_dict=edges_dict, 
        output_dir=satute_output_dir, 
        path_iqtree=path_iqtree, 
        alpha=alpha, 
        model=model)

    # Specify the path to output for the analysis
    output_dir = os.path.join(current_directory,"results_per_category_analysis", data_name)
    os.makedirs(output_dir, exist_ok=True)

    # Using SatuTe output for the analysis
    for edge_name in edges_dict.keys():
        print(f"Generate directory for {edge_name.upper()}")

        # Create a directory for the gene
        edge_directory = os.path.join(satute_output_dir, edge_name)
        os.makedirs(edge_directory, exist_ok=True)

        summary_category_zscore_per_branch(edge_directory, results_dir=output_dir, data_name=edge_name)

   