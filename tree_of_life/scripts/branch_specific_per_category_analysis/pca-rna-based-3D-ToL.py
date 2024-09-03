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
    per_category_analysis,
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
        "branch_to_eukaryota": "(Node1612*,Node77*)",
        "branch_to_yeast": "(Eukaryota_Opisthokonta_Nucletmycea_Fungi_Dikarya_Ascomycota_Saccharomycotina_Saccharomycetes_Saccharomycetales_Saccharomyces_cerevisiae_baker_s_yeast,Node1652*)",
    }

    # Specify the path to your considered data 
    data_name = "rRNA_based_3D_tree"
    input_dir = os.path.join(current_directory, "data", data_name)
    fasta_file = find_file_with_suffix(data_name, ".fasta", input_dir)
    model="GTR+G4"
    

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
    per_category_analysis(input_dir, results_dir=output_dir)
