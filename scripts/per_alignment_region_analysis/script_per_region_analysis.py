import pandas as pd
import sys
import os

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils.script_handle_satute_components import (
    calculate_variance, 
    summarize_component_data,
)

from utils.script_analyses_utils import (
    get_regions_from_annotation
)


# def calculate_zscore_per_region(branch_df, variance, region_info):

#     zscore_dict = {}
#     data = branch_df['coherence'].values

#     for region, sites  in region_info.items():
#         nsites = len(sites)
#         if nsites == 0:  # Handle cases where a region might have no sites
#             zscore_dict[region] = np.nan
#             continue

#         average  = sum(data[sites])/nsites
#         zscore_dict[region]= average/np.sqrt(variance/len(sites))

#     return zscore_dict




def calculate_zscore_per_region(branch_df, variance, region_info):
    zscore_dict = {}

    # Create a mapping from site values to their indices in the DataFrame
    site_to_index = {site: idx for idx, site in enumerate(branch_df['site'].values)}
    
    # Extract coherence values as a NumPy array
    data = branch_df['coherence'].values  

    for region, sites in region_info.items():
        # Get the actual indices in the DataFrame for the sites in this region
        site_indices = [site_to_index[site] for site in sites if site in site_to_index]
        nsites = len(site_indices)
        
        if nsites == 0:  # Handle cases where a region might have no valid sites
            zscore_dict[region] = np.nan
            continue
        
        # Calculate the average coherence for the region
        average = np.mean(data[site_indices])
        
        # Calculate the Z-score
        zscore = average / np.sqrt(variance / nsites)
        zscore_dict[region] = zscore

    return zscore_dict



def calculate_region_zscores_per_branch( component_dataframe, region_info): 
    unique_branches = component_dataframe['branch'].unique()
    if len(unique_branches) != 1:
        print(f"Error: Data frame includes data for {len(unique_branches)} branches!")
        return pd.DataFrame()

    variance = calculate_variance(component_dataframe)
    region_score = calculate_zscore_per_region(component_dataframe, variance,region_info)
        
    return region_score, variance

def per_region_analysis(annotation_file, satute_input_dir, edge_list=None, results_dir=None, data_name =None):
    # If results_dir is None, use satute_input_dir
    if results_dir is None:
        results_dir = satute_input_dir
    if data_name is None: 
        data_name = os.path.basename(satute_input_dir)

    # Get regions information for annotation file
    region_info = get_regions_from_annotation(annotation_file)

    ### Summarize all data of the coeherence coefficients per site in the directory
    summarized_component_data_dict = summarize_component_data(satute_input_dir, results_dir, data_name)
    region_dict = {}
    global_variance_list = []

    for dataset in summarized_component_data_dict.keys():
        unique_branches = summarized_component_data_dict[dataset]['branch'].unique()
        branch_results_dict = {}  # Dictionary to hold rolling results for each branch

        for branch in unique_branches:
            if edge_list is None or branch in edge_list:
                # Filter the dataset for the current branch
                branch_data = summarized_component_data_dict[dataset][summarized_component_data_dict[dataset]['branch'] == branch]

                # Calculate custom rolling metric using the window size
                region_zscores, variance = calculate_region_zscores_per_branch(branch_data, region_info)
          
                if region_zscores:  # Only proceed if region_zscores is not empty
                    # Convert the region_zscores dictionary to a DataFrame
                    region_df = pd.DataFrame.from_dict(region_zscores, orient='index', columns=[branch])

                    # Reset the index to convert the region names from the index to a column
                    region_df = region_df.reset_index()

                    # Rename the index column to 'region'
                    region_df = region_df.rename(columns={'index': 'region'})
                    
                    # Store the DataFrame in branch_results_dict under the branch name
                    branch_results_dict[branch] = region_df 
                    global_variance_list.append({'dataset': dataset, 'branch': branch, 'variance': variance})
    print(branch_results_dict)

    if branch_results_dict:
        # Combine all branch DataFrames into a single DataFrame
        combined_region_df = pd.concat(branch_results_dict.values(), axis=1)

        # Ensure the 'region' column is not duplicated during concatenation
        combined_region_df = combined_region_df.loc[:,~combined_region_df.columns.duplicated()]
        
        # Store the combined DataFrame in region_dict under the dataset key
        region_dict[dataset] = combined_region_df

        # Define the path to save the results
        region_csv_path = os.path.join(results_dir, f"{dataset}_per_region_zscores.csv")

        # Save the combined DataFrame to a CSV file
        combined_region_df.to_csv(region_csv_path, index=False)

        plot_zscores_per_region(combined_region_df, dataset, results_dir) 
    else: 
        print("ohje")

    if global_variance_list:  # Check if global_variance_list is not empty
        # Convert the list of dictionaries into a DataFrame
        global_variance_df = pd.DataFrame(global_variance_list)
        
        # Save the aggregated variance data
        csv_file_path = os.path.join(results_dir, f"{data_name}_global_variance.csv")
        global_variance_df.to_csv(csv_file_path, index=False)

    return combined_region_df


def plot_zscores_per_region(region_df, dataset_name, results_dir, edge_list=None):
    # Calculate the maximum and minimum z-scores across all branches
    y_max = region_df.drop(columns=['region']).max().max()
    y_min = region_df.drop(columns=['region']).min().min()

    # If edge_list is not provided, include all branches (all columns except 'region_name')
    if edge_list is None:
        edge_list = region_df.columns.difference(['region']).tolist()

    output_pdf = os.path.join(results_dir, f"{dataset_name}_per_region_zscores.pdf")

    with PdfPages(output_pdf) as pdf:
        # Plot each branch on the same figure
        for branch in edge_list:
            # Filter the data for the current branch
            branch_data = region_df[['region', branch]].copy()
            branch_data = branch_data.rename(columns={branch: 'z-score'})

            # Check if branch_data is empty
            if branch_data['z-score'].isna().all():
                print(f"No data found for branch {branch}. Skipping...")
                continue

            # Create the plot
            plt.figure(figsize=(10, 6))

            # Bar plot for the z-scores
            sns.barplot(
                data=branch_data, 
                x='region', 
                y='z-score', 
                color='steelblue'
            )
            
            # Final plot adjustments
            plt.title(f"Branch: {branch}", size=16)
            plt.ylim(y_min - 0.15, y_max + 0.15)
            plt.grid(False)
            plt.tight_layout()

            # Save the current figure to the PDF
            pdf.savefig()
            plt.close()  # Close the figure to free up memory

    print(f"All plots have been saved to {output_pdf}")






if __name__ == "__main__":

    # Get the current working directory
    current_directory = os.getcwd()

    region_annotation_file = os.path.join(current_directory, "../../example/data/region_annotation.csv")
   
    
    """ Analysis for all branches """

    data_name = "example_all_branches"

    # Specify the path to your considered  input data 
    input_dir = os.path.join(current_directory, "../../example/SatuTe_with_rate_heterogeneity")


    # Specify the path to output
    output_dir = os.path.join(current_directory,"../../example/results_per_region_analysis/", data_name)
    os.makedirs(output_dir, exist_ok=True)

    per_region_analysis(region_annotation_file, satute_input_dir=input_dir, results_dir=output_dir)



    """ Analysis for specific branches """

    data_name = "example_specific_branches"

    # Specify the path to your considered  input data 
    input_dir = os.path.join(current_directory, "../../example/SatuTe_without_rate_heterogeneity")

    # Specify the path to output
    output_dir = os.path.join(current_directory,"../../example/results_per_region_analysis/", data_name)
    os.makedirs(output_dir, exist_ok=True)

    per_region_analysis(region_annotation_file, input_dir, edge_list=["(A1, Node1*)","(Node4*, Node2*)"], results_dir=output_dir)
