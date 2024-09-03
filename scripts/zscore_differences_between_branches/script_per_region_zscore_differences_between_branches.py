import pandas as pd
import itertools
import sys
import os

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from per_alignment_region_analysis.script_per_region_analysis import (
    per_region_analysis,
)

from utils.script_handle_data import(
    subdictonary_additional_info
)

def calculate_all_pairwise_differences(region_zscores_df):
    # Initialize a list to store the results
    results = []

    # Extract the list of branch names (columns) excluding the 'region_name' column
    branches = region_zscores_df.columns.difference(['region'])

    # Iterate over each region
    for index, row in region_zscores_df.iterrows():
        region = row['region']
        
        # Generate all pairwise combinations of branches
        for branch_one, branch_two in itertools.combinations(branches, 2):
            # Calculate the difference in z-scores between the two branches
            zscore_difference = row[branch_one] - row[branch_two]
            
            # Append the result to the list
            results.append({
                'region': region,
                'zscore_difference': zscore_difference,
                'branch_one': branch_one,
                'branch_two': branch_two
            })

    # Convert the list of results into a DataFrame
    result_df = pd.DataFrame(results)

    return result_df


def per_region_zscore_differences_branches(annotation_file, satute_input_dir, edge_list=None, results_dir=None, data_name =None):
    # If results_dir is None, use satute_input_dir
    if results_dir is None:
        results_dir = satute_input_dir
    if data_name is None: 
        data_name = os.path.basename(satute_input_dir)

    region_zscores_df = per_region_analysis(annotation_file, satute_input_dir, edge_list, results_dir, data_name)

    subdictonary_additional_info(results_dir)
    
    differences_df = calculate_all_pairwise_differences(region_zscores_df)
        
    csv_path = os.path.join(results_dir, f"{data_name}_per_region_zscore_differences.csv")
    differences_df.to_csv(csv_path, index=False)

    plot_zscore_differences_per_branch_pair(differences_df, data_name, results_dir)


def plot_zscore_differences_per_branch_pair(differences_df, dataset_name, results_dir, title=None):
    y_max = differences_df['zscore_difference'].max()
    y_min = differences_df['zscore_difference'].min()
    # Sort the differences in ascending order
    differences_df = differences_df.sort_values(by='zscore_difference')

    # Default to all unique branch pairs if edge_list is not provided
    edge_list = differences_df[['branch_one', 'branch_two']].drop_duplicates().values.tolist()

    if not title:
        output_pdf = os.path.join(results_dir, f"{dataset_name}_zscore_differences_per_branch_pair.pdf")
    else:
        output_pdf = os.path.join(results_dir, f"{title}_zscore_differences.pdf")

    with PdfPages(output_pdf) as pdf:
        # Plot each branch pair
        for branch_one, branch_two in edge_list:

            # Filter the data for the current branch pair
            branch_pair_data = differences_df[
                (differences_df['branch_one'] == branch_one) & (differences_df['branch_two'] == branch_two)
            ]

            # Check if branch_pair_data is empty
            if branch_pair_data.empty:
                print(f"No data found for branch pair {branch_one} - {branch_two}. Skipping...")
                continue

            # Create the plot
            plt.figure(figsize=(10, 6))

            # Bar plot for the z-score differences
            sns.barplot(
                data=branch_pair_data, 
                x='region', 
                y='zscore_difference', 
                palette=['steelblue' if x < 0 else 'red' for x in branch_pair_data['zscore_difference']],
            )
            
            # Final plot adjustments
            if not title: 
                title_text = f"Z-score Differences: {branch_one} vs {branch_two}"[:50]  # Limit title to 50 characters
            else: 
                title_text =  f"Z-score Differences: {title}"
            plt.title(title_text, size=16)
            plt.ylim(y_min - 0.15, y_max + 0.15)
            plt.grid(False)
            plt.tight_layout()

            plt.xlabel("Region", size=14)
            plt.ylabel("Z-score differences", size=14)

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
    output_dir = os.path.join(current_directory,"../../example/results_zscore_differences_branches/", data_name)
    os.makedirs(output_dir, exist_ok=True)

    per_region_zscore_differences_branches(region_annotation_file, satute_input_dir=input_dir, results_dir=output_dir)



    """ Analysis for specific branches """

    data_name = "example_specific_branches"

    # Specify the path to your considered  input data 
    input_dir = os.path.join(current_directory, "../../example/SatuTe_without_rate_heterogeneity")

    # Specify the path to output
    output_dir = os.path.join(current_directory,"../../example/results_zscore_differences_branches/", data_name)
    os.makedirs(output_dir, exist_ok=True)

    per_region_zscore_differences_branches(region_annotation_file, input_dir, edge_list=["(A1, Node1*)","(Node4*, Node2*)"], results_dir=output_dir)
