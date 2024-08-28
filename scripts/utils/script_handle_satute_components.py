import pandas as pd
import sys
import os

# Add the root directory (scripts) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from utils.script_handle_data import (
    process_gene_directory_components
)


def summarize_component_data(results_dir, output_dir, data_name=None):
    if data_name is None: 
        data_name = os.path.basename(results_dir)

    summarized_component_data = {}

    print("Summarize coherence coefficient data:")
    
    # List all subdirectories in the results_dir
    subdirs = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]

    if subdirs:  # If there are subdirectories, process each one
        for dataset_name in subdirs:
            print(f"Processing subfolder: {dataset_name}")
            gene_directory = os.path.join(results_dir, dataset_name)
            
            # Analysis summary for specific gene directory
            gene_component_data = process_gene_directory_components(gene_directory, dataset_name)
            
            if not gene_component_data.empty:
                # Sort the data by the "site" column and append to the summarized data
                gene_component_data.sort_values(by=["branch", "site"], inplace=True)
                summarized_component_data[dataset_name] = gene_component_data

    
    # If there are no subdirectories, process the results_dir itself
    print(f"Processing: {results_dir}")
    gene_component_data = process_gene_directory_components(results_dir, os.path.basename(results_dir))
    
    if not gene_component_data.empty:
        gene_component_data.sort_values(by=["branch", "site"], inplace=True)
        summarized_component_data[os.path.basename(results_dir)] = gene_component_data

    if summarized_component_data:  # Checking if the dictionary is not empty
        # Aggregate all DataFrames into a single DataFrame
        all_components_data = pd.concat(summarized_component_data.values(), ignore_index=True)
        
        # Save the aggregated data
        csv_file_path = os.path.join(output_dir, f"{data_name}_summarized_component_data.csv")
        all_components_data.to_csv(csv_file_path, index=False)

    print("")
    return summarized_component_data


def calculate_variance(component_dataframe):
    variance = 0
    rate_counts = component_dataframe['rate_category'].value_counts()
    sequence_len = sum(rate_counts)
    variance_by_rate = component_dataframe.groupby('rate_category')['category_variance'].first()
    for rate, count in rate_counts.items():
        variance += variance_by_rate[rate] * count/sequence_len
    return variance