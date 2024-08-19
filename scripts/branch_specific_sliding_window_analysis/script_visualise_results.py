import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import sys
import os

# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils.script_analyses_utils import (
    find_file_with_suffix_in_directory,
)

def plot_coherence_distribution_per_branch(results_dir, start_site=None, end_site=None, edge_list=None, bins = 10):
    """
    Plots the coherence coefficients of all considered site (coherence distribution) for each branch.

    Args:
        results_dir (str): Directory containing the result CSV file.
        start_site (int, optional): Start of the site range. If None, starts from the first site.
        end_site (int, optional): End of the site range. If None, ends at the last site.
        edge_list (list, optional): List of branches to include in the plot. If None, include all branches.
        bins (int, optional): Number of bins to use in the histogram. Default is 10.
    """
    # Load the result DataFrame
    component_dataframe = find_file_with_suffix_in_directory("_summarized_component_data.csv", results_dir)
    result_df= pd.read_csv(component_dataframe)

    # Filter by edge_list if provided
    if edge_list is not None:
        result_df = result_df[result_df['branch'].isin(edge_list)]

    output_pdf = os.path.join(results_dir, f"Coherence_coeffiecient_histogram_per_branch_{bins}_bins.pdf")
    # Create a PdfPages object to save multiple plots in a single PDF
    with PdfPages(output_pdf) as pdf:
        # Plot histogram for Coefficient column for each branch
        branches = result_df['branch'].unique()
        for branch in branches:
            branch_df = result_df[result_df['branch'] == branch]
            
            plt.figure()  # Create a new figure for each branch
            plt.hist(branch_df['coherence'], bins=bins, alpha=0.7, edgecolor='black')
            plt.xlabel('Coherence Coefficients', fontsize=12)
            plt.ylabel('Frequency', fontsize=12)
            plt.title(f'Coherence Histogram for Branch {branch}', fontsize=14)
            plt.grid(True)
            
            # Save the current figure to the PDF
            pdf.savefig()
            plt.close()  # Close the figure to free up memory

    print(f"All plots have been saved to {output_pdf}")



def plot_sliding_window_analysis_combined(csv_file, results_dir, edge_list=None):
    """
    Plots the sliding window analysis for each branch in the provided results file on the same plot and saves it in a PDF.

    Args:
        csv_file (str): Path to the CSV file containing the sliding window analysis data.
        results_dir (str): Directory where the output PDF will be saved.
        edge_list (list, optional): List of branches to include in the plot. If None, include all branches.
    """
    # Load the result DataFrame
    result_df = pd.read_csv(csv_file)

    # Filter by edge_list if provided
    if edge_list is not None:
        result_df = result_df[edge_list]
    else:
        edge_list = result_df.columns.tolist()

    output_pdf = os.path.join(results_dir, f"Sliding_window_analysis_combined_for_{len(edge_list)}_branches.pdf")

    # Create a new figure for the combined plot
    plt.figure(figsize=(10, 6))

    # Plot each branch on the same figure
    for branch in edge_list:
        plt.plot(result_df[branch], label=branch)
    
    # Add labels, title, and legend
    plt.xlabel('Site', fontsize=12)
    plt.ylabel('Z-score', fontsize=12)
    plt.title(f'Sliding Window Analysis Combined for {len(edge_list)} Branches', fontsize=14)
    plt.legend(title='Branch', loc='best')
    plt.grid(True)
    
    # Save the figure to a PDF
    plt.savefig(output_pdf)
    plt.close()  # Close the figure to free up memory

    print(f"Combined plot has been saved to {output_pdf}")






if __name__ == "__main__":

    # Get the current working directory
    current_directory = os.getcwd()

    # Specify the path to results
    results_dir = os.path.join(current_directory,"../../example/results_sliding_window_analysis/example_specific_branches")
    
    """" Visualise Coherence Coeffiecient Distribution per Branch """

    plot_coherence_distribution_per_branch(results_dir, edge_list=["(A1, Node1*)","(Node4*, Node2*)"], bins=20)


    """" Visualise Sliding window analysis for different Branches"""
    
    csv_file = find_file_with_suffix_in_directory("sliding_window_size_36.csv", results_dir)
    plot_sliding_window_analysis_combined(csv_file, results_dir)
   


   

