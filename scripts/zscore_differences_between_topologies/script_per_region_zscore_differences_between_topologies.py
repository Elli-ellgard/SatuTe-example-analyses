import pandas as pd
import sys
import os

# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from per_alignment_region_analysis.script_per_region_analysis import (
    per_region_analysis,
)

from zscore_differences_between_branches.script_per_region_zscore_differences_between_branches import (
    calculate_all_pairwise_differences,
    plot_zscore_differences_per_branch_pair
)

from utils.script_handle_data import(
    subdictonary_additional_info,
)
 
def combine_dict_to_df(data_dict):
    # Initialize an empty DataFrame for the final combined result
    combined_df = None

    # Loop through the dictionary to rename columns and merge the DataFrames
    for key, df in data_dict.items():
        # Rename columns
        renamed_df = df.rename(columns={col: f'{key}_{col}' for col in df.columns if col != 'region'})
        
        # Merge DataFrames on 'region', starting with the first DataFrame
        if combined_df is None:
            combined_df = renamed_df
        else:
            combined_df = pd.merge(combined_df, renamed_df, on='region')

    # Display the combined DataFrame
    return combined_df

def per_region_zscores_toplogies(region_zscores_dict, results_dir):
    subdictonary_additional_info(results_dir)
    zscore_dataframe = combine_dict_to_df(region_zscores_dict) 
    differences_df = calculate_all_pairwise_differences(zscore_dataframe)
        
    csv_path = os.path.join(results_dir, f"Per_region_zscore_differences_topologies.csv")
    differences_df.to_csv(csv_path, index=False)

    plot_zscore_differences_per_branch_pair(differences_df, "Different_topologies", results_dir)


if __name__ == "__main__":

    # Get the current working directory
    current_directory = os.getcwd()

    # Specify the path to output
    output_dir = os.path.join(current_directory,"../../example/results_zscore_differences_topologies/")
    os.makedirs(output_dir, exist_ok=True)

    region_annotation_file = os.path.join(current_directory, "../../example/data/region_annotation.csv")
    
    region_zscores = {}

    # original tree analysis / input files
    data_name = "true_tree"
    satute_input_dir = os.path.join(current_directory, "../../example/SatuTe_without_rate_heterogeneity")
    edges_dict= {
        "branch": "(Node6*, Node5*)",
    }
    region_zscores[data_name] = per_region_analysis(region_annotation_file, satute_input_dir,  list(edges_dict.values()), output_dir, data_name)

    # random tree analysis
    data_name = "random_tree"
    satute_input_dir = os.path.join(current_directory, "../../example/comparison_topology")
    edges_dict= {
        "branch": "(Node5*, Node3*)",
    }
    region_zscores[data_name] = per_region_analysis(region_annotation_file, satute_input_dir,  list(edges_dict.values()), output_dir, data_name)

    per_region_zscores_toplogies(region_zscores, output_dir)

    












