import pandas as pd

import sys
import os

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils.script_handle_data import (
    extract_relative_category_rates_from_satute
)


# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils.script_handle_data import (
    process_directory,
    summarize_z_scores_categories,
)

from utils.script_analyses_utils import (
    find_file_with_suffix_in_directory,
)

from utils.script_handle_tree import (
    get_newick_string_from_satute,
)

def summary_analysis_branch_saturation(input_dir, data_name, newick_string,results_dir):

    if os.path.isdir(input_dir):
        print(f"Processing folder: {data_name}") 

        # Analysis summary for specific gene directory
        process_directory(input_dir, data_name, newick_string, results_dir)
    else: 
        print("Folder does not exists!")


def custom_transform(y, y_max):
    """
    Custom transformation function that scales the input values.

    Args:
        y (pd.Series or np.array): The input values to transform.
        y_max (float): The maximum value to scale the input values to.

    Returns:
        pd.Series or np.array: The transformed values.
    """
    y_transformed = (y / y.max()) * y_max
    return y_transformed

    
def plot_summary_z_scores_per_branch(z_score_data, category_rates, dataset_name, results_dir, edge_list=None):
    # Calculate the maximum and minimum z-scores
    y_max = max(z_score_data["z_score"])
    y_min = min(min(z_score_data["z_score"]),0)

    # If edge_list is not provided, include all branches
    if edge_list is None:
        edge_list = z_score_data["branch"].unique()

    # Create a dictionary to map rate categories to relative rates
    category_rate_map = dict(zip(category_rates["Category"], category_rates["Relative_rate"]))

    # Ensure the rate categories have a consistent order
    category_order = category_rates["Category"].tolist()
    z_score_data["rate_category"] = pd.Categorical(z_score_data["rate_category"], categories=category_order, ordered=True)

    output_pdf = f"{results_dir}/{dataset_name}_summary_zscores.pdf"

    with PdfPages(output_pdf) as pdf:
        # Plot each branch on the same figure
        for branch in edge_list:

            # Filter the data for the current branch
            branch_data = z_score_data[z_score_data["branch"] == branch]

                        # Check if branch_data is empty
            if branch_data.empty:
                print(f"No data found for branch {branch}. Skipping...")
                continue

            # Check if the columns 'z_alpha' and 'z_alpha_bonferroni_corrected' exist
            if 'z_alpha' not in branch_data.columns or 'z_alpha_bonferroni_corrected' not in branch_data.columns:
                print(f"Missing z_alpha data for branch {branch}. Skipping...")
                continue

            # Calculate z_alpha and z_alpha_corrected for this branch
            z_alpha =   branch_data["z_alpha"].iloc[0]
            z_alpha_corrected = branch_data["z_alpha_bonferroni_corrected"].iloc[0]

            # Create a copy of the DataFrame for transformation
            branch_data_transformed = branch_data.copy()
            branch_data_transformed['number_informative_sites_transformed'] = custom_transform(
                branch_data['number_of_sites'], y_max
            )

            # Create the plot
            plt.figure(figsize=(10, 6))

            # Bar plot for the transformed number of informative sites
            sns.barplot(
                data=branch_data_transformed, 
                x='rate_category', 
                y='number_informative_sites_transformed', 
                color='lightgrey',
                order=category_order
            )
            
            # Overlay points and lines for z-scores
            sns.pointplot(
                data=branch_data, 
                x='rate_category', 
                y='z_score', 
                hue='branch', 
                markers='o', 
                linestyles='dotted',
                order=category_order
            )
            plt.legend().remove()

            # Add horizontal lines for z_alpha and z_alpha_corrected
            plt.axhline(y=z_alpha, color='black', linestyle='dashed')
            plt.text(len(category_rates) - 0.1, z_alpha + 1, r'$z_{\alpha}$', verticalalignment='center', size=12)

            plt.axhline(y=z_alpha_corrected, color='red', linestyle='dashed')
            plt.text(len(category_rates) - 0.1, z_alpha_corrected + 1.2, r'$z_{\alpha[adjust]}$', color='red', size=12)

            # Customize x-ticks with relative rates as labels
            relative_rate_labels = [f"{category_rate_map.get(lbl, lbl):.4f}" for lbl in category_order]
            plt.xticks(ticks=range(len(relative_rate_labels)), labels=relative_rate_labels, size=10)

            # Customize axis labels and limits
            plt.xlabel("Rate category", size=14)
            plt.ylabel("Z-score", size=14)
            plt.ylim(y_min -1, y_max +1)

            # Secondary y-axis for sequence length
            sec_ax = plt.gca().twinx()
                
            sec_ax.set_ylim((y_min-1)*branch_data["number_of_sites"].max()/(y_max -1), 
                    branch_data_transformed['number_informative_sites_transformed'].max()/(y_max-1)*branch_data["number_of_sites"].max())
            
            sec_ax.set_ylabel("Sequence length", size=14, color='darkgrey')
            sec_ax.tick_params(axis='y', labelsize=10, colors='darkgrey')

            # Final plot adjustments
            title_text = f"Branch: {branch}"[:50]  # Limit title to 50 characters
            plt.title(title_text, size=16)
            plt.grid(False)
            plt.tight_layout()

            # Save the current figure to the PDF
            pdf.savefig()
            plt.close()  # Close the figure to free up memory

    print(f"All plots have been saved to {output_pdf}")


def summary_category_zscore_per_branch(satute_input_dir, data_name, results_dir, edge_list=None):

    zscore_data = summarize_z_scores_categories(satute_input_dir, data_name)

    satute_file = find_file_with_suffix_in_directory(".satute", satute_input_dir)
    category_rates = extract_relative_category_rates_from_satute(satute_file)
    
    plot_summary_z_scores_per_branch(zscore_data, category_rates, data_name, results_dir, edge_list)
 


def per_category_analysis(satute_input_dir, edge_list=None, results_dir=None, data_name =None):
    # If results_dir is None, use satute_input_dir
    if results_dir is None:
        results_dir = satute_input_dir
    if data_name is None: 
        data_name = os.path.basename(satute_input_dir)

    #get considered tree 
    print("Summary Branch Saturation\n")
    satute_file = find_file_with_suffix_in_directory(".satute", satute_input_dir)
    newick_string = get_newick_string_from_satute(satute_file)
    summary_analysis_branch_saturation(satute_input_dir, data_name, newick_string, results_dir)

    print("\nSummary Zscores")
    summary_category_zscore_per_branch(satute_input_dir, data_name, results_dir, edge_list)







if __name__ == "__main__":

    # Get the current working directory
    current_directory = os.getcwd()

    # Specify the path to your considered  input data 
    input_dir = os.path.join(current_directory, "../../example/SatuTe_with_rate_heterogeneity")

    
    """ Analysis for all branches """

    data_name = "example_all_branches"

    # Specify the path to output
    output_dir = os.path.join(current_directory,"../../example/results_per_category_analysis/", data_name)
    os.makedirs(output_dir, exist_ok=True)

    per_category_analysis(satute_input_dir=input_dir, results_dir=output_dir)



    """ Analysis for specific branches """

    data_name = "example_specific_branches"

    # Specify the path to output
    output_dir = os.path.join(current_directory,"../../example/results_per_category_analysis/", data_name)
    os.makedirs(output_dir, exist_ok=True)

    per_category_analysis(input_dir, edge_list=["(A1, Node1*)","(Node4*, Node2*)"], results_dir=output_dir)

