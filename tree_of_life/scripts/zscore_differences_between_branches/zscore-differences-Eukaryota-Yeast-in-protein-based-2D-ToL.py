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
    model = "LG+G4"

    folder_dict = {}

    # original tree analysis / input files
    dataname = "protein_based_2D_tree"
    input_dir = os.path.join(current_directory, "data", dataname)
    region_annotation_file = os.path.join(input_dir, "protein_annotation.csv")
    edges_dict= {
        "branch_to_eukaryota": "(Node2691*, Node2692*)",
        "branch_to_yeast": "(Eukaryota_Opisthokonta_Fungi_Dikarya_Ascomycota_saccharomyceta_Saccharomycotina_Saccharomycetes_Saccharomycetales_Saccharomycetaceae_Saccharomyces_cerevisiae, Node2721*)",
    }

    # Preprocessing: generate directories + run SatuTe and IQ-Tree
    satute_dir= os.path.join(current_directory,"SatuTe_results", dataname)
    os.makedirs(satute_dir, exist_ok=True) 

    branch_specific_preprocessing(
        input_dir=input_dir, 
        edges_dict=edges_dict, 
        output_dir=satute_dir, 
        path_iqtree=path_iqtree, 
        alpha=alpha, 
        model =model)

    # summarize component data for the considered edges
    folder_dict = {}
    for edge in edges_dict.keys():
        folder_dict[edge] = os.path.join(satute_dir, edge) 
  
   
    ''' Comparison region zscores '''

    # Specify the path to output
    output_dir = os.path.join(current_directory,"results_zscore_differences_branches/", dataname)
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

    per_region_zscores_toplogies(region_zscores, output_dir, title="protein_based_2D_Eukaryota_minus_Yeast")

