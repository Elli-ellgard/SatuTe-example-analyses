import pandas as pd
import numpy as np
import sys
import os

# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from utils.script_handle_satute_components import (
    calculate_variance, 
    summarize_component_data,
)

def calculate_window_zscore(window_data,variance):
    window_size = len(window_data)
    weights = np.linspace(1, 0.5, window_size)  # Linearly decreasing weights
    weighted_average = np.average(window_data, weights=weights)
    return weighted_average/np.sqrt(variance/window_size)


def calculate_window_zscores_per_branch( component_dataframe, window_size): 

    unique_branches = component_dataframe['branch'].unique()
    if len(unique_branches) != 1:
        print(f"Error: Data frame includes data for {len(unique_branches)} branches!")
        return pd.DataFrame()

    data = np.array(component_dataframe['coherence'])
    variance = calculate_variance(component_dataframe)
    window_score = pd.Series(data).rolling(window=window_size).apply(lambda x: calculate_window_zscore(x, variance), raw=True)
    
    # Shift the results to align them to the middle of the window
    middle_position_shift = -(window_size // 2)
    window_score_centered = window_score.shift(middle_position_shift)
    
    return window_score_centered, variance
         

def sliding_window_analysis(satute_input_dir, window_size, edge_list=None, results_dir=None, data_name =None):
    # If results_dir is None, use satute_input_dir
    if results_dir is None:
        results_dir = satute_input_dir
    if data_name is None: 
        data_name = os.path.basename(satute_input_dir)

    ### Summarize all data of the coeherence coefficients per site in the directory
    summarized_component_data_dict = summarize_component_data(satute_input_dir, results_dir, data_name)
    sliding_window_dict = {}
    global_variance_list = []

    for dataset in summarized_component_data_dict.keys():
        unique_branches = summarized_component_data_dict[dataset]['branch'].unique()
        branch_results_dict = {}  # Dictionary to hold rolling results for each branch

        for branch in unique_branches:
            if edge_list is None or branch in edge_list:
                # Filter the dataset for the current branch
                branch_data = summarized_component_data_dict[dataset][summarized_component_data_dict[dataset]['branch'] == branch]

                # Calculate custom rolling metric using the window size
                window_score_centered, variance = calculate_window_zscores_per_branch(branch_data, window_size)

                if not window_score_centered.empty:
                    # Store the rolling centered data in the dictionary under the branch name
                    branch_results_dict[branch] = window_score_centered
                    # Also store the variance for that branch
                    global_variance_list.append(
                        {'dataset': dataset, 
                         'branch': branch, 
                         'global_variance': variance,
                        })

        if branch_results_dict:
            # Convert the dictionary of rolling results into a DataFrame
            sliding_window_df = pd.DataFrame(branch_results_dict)
            sliding_window_dict[dataset] = sliding_window_df

            # Save the sliding window results for the current dataset
            sliding_window_csv_path = os.path.join(results_dir, f"{dataset}_sliding_window_size_{window_size}.csv")
            sliding_window_df.to_csv(sliding_window_csv_path, index=False)

    if global_variance_list:  # Check if global_variance_list is not empty
        # Convert the list of dictionaries into a DataFrame
        global_variance_df = pd.DataFrame(global_variance_list)
        
        # Save the aggregated variance data
        csv_file_path = os.path.join(results_dir, f"{data_name}_global_variance.csv")
        global_variance_df.to_csv(csv_file_path, index=False)


if __name__ == "__main__":

    # Get the current working directory
    current_directory = os.getcwd()


    """ Sliding window analysis for all branches """

    data_name = "example_all_branches"

    # Specify the path to your considered  input data 
    input_dir = os.path.join(current_directory, "../../example/SatuTe_with_rate_heterogeneity")

    # Specify the path to output
    output_dir = os.path.join(current_directory,"../../example/results_sliding_window_analysis/", data_name)
    os.makedirs(output_dir, exist_ok=True)

    sliding_window_analysis(satute_input_dir=input_dir, window_size=40, results_dir=output_dir)



    """ Sliding window analysis for specific branches """

    data_name = "example_specific_branches"

    # Specify the path to your considered  input data 
    input_dir = os.path.join(current_directory, "../../example/SatuTe_without_rate_heterogeneity")

    # Specify the path to output
    output_dir = os.path.join(current_directory,"../../example/results_sliding_window_analysis/", data_name)
    os.makedirs(output_dir, exist_ok=True)

    sliding_window_analysis(input_dir, window_size=36, edge_list=["(A1, Node1*)","(Node4*, Node2*)"], results_dir=output_dir)




