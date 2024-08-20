import os
import sys
# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from utils.script_analyses_utils import (
    find_file_with_suffix, 
    find_files_with_suffix_in_directory
)

from utils.script_preprocessing_and_generating_satute_output import (
    branch_specific_preprocessing,
)

from branch_specific_sliding_window_analysis.script_sliding_window_analysis import (
    sliding_window_analysis,
)

from branch_specific_sliding_window_analysis.script_visualise_results import (
    plot_sliding_window_analysis_combined,
)

   

if __name__ == "__main__":

    # Set paths to IQ-TREE and Python executable
    path_iqtree = "iqtree2"

    # Set significance level
    alpha = str(0.01)

    """ Set directories and input files """
    # Get the current working directory
    current_directory = os.path.join(os.getcwd(), "../../", "tree_of_life/")


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
    

    # Specify the path to output filtered fasta file
    output_dir = os.path.join(current_directory,"results_sliding_window_analysis", data_name)
    os.makedirs(output_dir, exist_ok=True)

    # Preprocessing: generate directories + run SatuTe and IQ-Tree
    branch_specific_preprocessing(input_dir, edges_dict, output_dir, path_iqtree, alpha, model)

    # Using SatuTe output for the analysis
    sliding_window_analysis(output_dir, window_size=36)

    csv_files = find_files_with_suffix_in_directory("sliding_window_size_36.csv", output_dir)
    plot_sliding_window_analysis_combined(csv_files, output_dir)
    