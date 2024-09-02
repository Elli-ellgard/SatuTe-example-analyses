import os
import re
import pandas as pd
import sys
import io
from pandas import DataFrame
import shutil
from utils.script_handle_tree import  write_nexus_file


def subdictonary_additional_info(current_dir):
    directory_name = os.path.join(current_dir, "additional_info")

    # Create the directory if it doesn't exist
    os.makedirs(directory_name, exist_ok=True)

    # List of file extensions to move
    file_extensions = ['.csv', '.pdf']

    # Move files with the specified extensions
    for file_name in os.listdir(current_dir):
        if any(file_name.endswith(ext) for ext in file_extensions):
            shutil.move(os.path.join(current_dir, file_name), os.path.join(directory_name, file_name))

   

def extract_relative_category_rates_from_satute(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    string_table = []
    # Find the start of the table
    table = False
    for  line in lines:
        if "Category" in line and "Relative_rate" in line:         
            table = True
            continue
            
        if table: 
            if "SPECTRAL DECOMPOSITION" in line: 
                break
            if line.strip()!= "":
                string_table.append(line.strip())

    if not table:
        raise ValueError("Table not found in the file.")
    
    # Join the lines into a single string and convert to a DataFrame
    table_str = '\n'.join(string_table)
    columns = [
        "Category", 
        "Relative_rate", 
        "Proportion", 
        "Empirical_Proportion"
    ]
    df = pd.read_csv(io.StringIO(table_str), delim_whitespace=True, names=columns)    

    # Select only the first two columns
    df = df[["Category", "Relative_rate"]]

    # Rename the Category column to c1, c2, ...
    df["Category"] = ['c{}'.format(i + 1) for i in range(len(df))]
    
    return df

def parse_file_to_data_frame(file_path):
    try:
        # Read the file into a dataframe
        df = pd.read_csv(file_path, delimiter="\t")

        return df

    except FileNotFoundError:
        raise Exception(f"File not found: {file_path}")



def get_column_names_with_prefix(data_frame, prefix):
    # Filter the columns using the specified prefix
    columns_with_prefix = data_frame.columns[
        data_frame.columns.str.startswith(prefix)
    ].tolist()
    return columns_with_prefix

def build_categories_by_sub_tables(data_frame: DataFrame):
    rate_category_dictionary = {}

    # Assuming you already have a dataframe called 'dataframe'
    # Call the get_columns_with_prefix function to retrieve columns with a specific prefix
    prefix = "p"  # Specify the desired prefix

    columns_with_prefix = get_column_names_with_prefix(data_frame, prefix)

    # Create dictionaries using column names with the specified prefix as names
    rate_category_dictionary = {column: [] for column in columns_with_prefix}

    for index, row in data_frame.iterrows():
        p_row = row.filter(like="p")

        rate_category_dictionary[p_row.idxmax()].append(int(row["Site"]) - 1)

    return rate_category_dictionary

def get_sites_per_category(directory_path):
    file_suffix = ".siteprob"
    log_files = [f for f in os.listdir(directory_path) if f.endswith(file_suffix)]
    if log_files:
        log_file = log_files[0]
        site_probability =  parse_file_to_data_frame(os.path.join(directory_path, log_file))
        categorized_sites = build_categories_by_sub_tables(site_probability) 
        return(categorized_sites)
    else: 
        return(-1)




def extract_file_name_without_suffix(file_path):
    # Use os.path.basename to get the file name from the path
    file_name = os.path.basename(file_path)
    
    # Use os.path.splitext to split the file name and extension
    file_name_without_suffix, _ = os.path.splitext(file_name)

    # Loop to remove multiple suffixes
    while True:
        file_name_without_suffix, suffix = os.path.splitext(file_name_without_suffix)
        if not suffix:
            break
    return file_name_without_suffix

def extract_rows_with_value(df, column_name, target_value):
    filtered_rows = df[df[column_name] == target_value]
    return filtered_rows

def rearrange_dataframe(summarized_data):

    # Rearrange the data frame
    rearranged_df = summarized_data.groupby('branch').agg(
        category=pd.NamedAgg(column='rate_category', aggfunc=lambda x: list(x.unique())),
        result_test=pd.NamedAgg(column='decision_test', aggfunc=lambda x: list(x.unique())),
        dataset=pd.NamedAgg(column='dataset', aggfunc=lambda x: list(x.unique())),
    ).reset_index()
   
    return  rearranged_df

def rearrange_dataframe_correction(summarized_data):

    # Rearrange the data frame
    rearranged_df = summarized_data.groupby('branch').agg(
        category=pd.NamedAgg(column='rate_category', aggfunc=lambda x: list(x.unique())),
        result_test=pd.NamedAgg(column='decision_bonferroni_corrected', aggfunc=lambda x: list(x.unique())),
        dataset=pd.NamedAgg(column='dataset', aggfunc=lambda x: list(x.unique())),
    ).reset_index()
   
    return  rearranged_df


def extract_rows_of_csv(data):
    # Replace this function with your specific action to extract rows with decision "Saturated"
    return data[data['decision_test'] == 'Saturated'].copy()

def extract_rows_of_csv_correction(data):
    # Replace this function with your specific action to extract rows with decision "Saturated"
    return data[data['decision_bonferroni_corrected'] == 'Saturated'].copy()

def add_seq_len_category_from_siteprob(data, categorized_sites, directory_path):
    category = list(data.loc[:, 'category_rate'])[0]
    seq_len = len(categorized_sites[category])
    new_value =f"{category}({seq_len})"
    data.loc[:, 'category_rate'] = new_value
    return data

def add_seq_len_category_from_log(data, directory_path):
    file_suffix = ".satute.log"
    log_file = [f for f in os.listdir(directory_path) if f.endswith(file_suffix)][0]
       
    with open(os.path.join(directory_path, log_file), 'r') as file:
        content = file.read()

    category = list(data.loc[:, 'category_rate'])[0]

    pattern = re.compile(fr'{re.escape(category)}, Site per category (\d+)')
    match = pattern.search(content)

    if match:
        seq_len = int(match.group(1))
        new_value =f"{category}({seq_len})"
        data.loc[:, 'category_rate'] = new_value
                
    else:
        print(f"No match found for {category}.")
    
    return data


def summarize_saturated_results_categories(directory_path, dataset_name):
    # Initialize an empty DataFrame to store the summarized data
    summarized_data = pd.DataFrame()
    summarized_data_corrrection = pd.DataFrame()

    # Get a list of satute results csv files
    file_suffix = ".satute.csv"
    satute_result_files = [f for f in os.listdir(directory_path) if f.endswith(file_suffix)]

    # Summarize the saturated results
    for file in satute_result_files:
        # Your specific action goes here
        print("Processing csv file:", file)
        
        # Read CSV file
        data = pd.read_csv(os.path.join(directory_path, file))
        
        # Extract rows with decision "Saturated"
        saturated_rows = extract_rows_of_csv(data)
        saturated_rows_correction = extract_rows_of_csv_correction(data)

        if not saturated_rows.empty:
            # Add dataset name
            #saturated_rows['dataset'] = dataset_name
            saturated_rows.loc[:, 'dataset'] = dataset_name
            saturated_rows_correction.loc[:, 'dataset'] = dataset_name


            # #add seq_length to category
            # categorized_sites = get_sites_per_category(directory_path)
            # saturated_rows= add_seq_len_category_from_siteprob(saturated_rows, categorized_sites, directory_path)
            
            # Append the filtered rows to the summarized data
            summarized_data = pd.concat([summarized_data, saturated_rows], ignore_index=True)
            summarized_data_corrrection = pd.concat([summarized_data_corrrection, saturated_rows_correction], ignore_index=True)

    return summarized_data, summarized_data_corrrection

def summarize_z_scores_categories(directory_path, dataset_name):
    # Initialize an empty DataFrame to store the summarized data
    zscore_data = pd.DataFrame()

    # Get a list of satute results csv files
    file_suffix = ".satute.csv"
    satute_result_files = [f for f in os.listdir(directory_path) if f.endswith(file_suffix)]

    # Summarize the saturated results
    for file in satute_result_files:
        # Your specific action goes here
        print("Processing csv file:", file)
        
        # Read CSV file
        data = pd.read_csv(os.path.join(directory_path, file))
        
        
        filtered_data = data[['branch', 'z_score', 'z_alpha','z_alpha_bonferroni_corrected','number_of_sites', 'rate_category']].copy()
        if not filtered_data.empty:
            # Add dataset name
            filtered_data.loc[:, 'dataset'] = dataset_name
            zscore_data = pd.concat([zscore_data, filtered_data])

    return zscore_data

def summarize_components_categories(directory_path, dataset_name):
    # Initialize an empty DataFrame to store the summarized data
    component_data = pd.DataFrame()

    # Get a list of satute results csv files
    file_suffix = ".components.csv"
    satute_result_files = [f for f in os.listdir(directory_path) if f.endswith(file_suffix)]

    # Summarize the saturated results
    for file in satute_result_files:
        # Your specific action goes here
        print("Processing csv file:", file)
        
        # Read CSV file
        data = pd.read_csv(os.path.join(directory_path, file))
        #filtered_data = data[data['test_statistic'] != 'duplicate sequence']
        
        filtered_data = data 
        if not filtered_data.empty:
            # Add dataset name
            filtered_data.loc[:, 'dataset'] = dataset_name
            component_data = pd.concat([component_data, filtered_data])

    return component_data

def process_gene_directory(directory_path, dataset_name,newick_string):

    # Data frame to store the summarized results for saturation over all categories
    summarized_data, summarized_data_correction = summarize_results_categories(directory_path, dataset_name)
    
    if not summarized_data.empty:  
        # Rearrange summarized data
        rearranged_data = rearrange_dataframe(summarized_data)

        # Save rearranged data
        csv_file_path = os.path.join(directory_path, f"{dataset_name}_summarized_data.csv")
        # Write DataFrame to CSV file
        rearranged_data.to_csv(csv_file_path, index=False)

        #Map rearranged data to newick string 
        if newick_string:
            # map rearranged data to it
            #print(newick_string)
            result_nexus_file = os.path.join(directory_path, f"{dataset_name}_summary_tree.nex")
            write_nexus_file(newick_string,result_nexus_file,rearranged_data)
    else:
        print("No saturated branches for ", dataset_name, "!")
       
    return summarized_data

def process_directory(directory_path, dataset_name,newick_string, results_dir= None):
    if results_dir is None:
        results_dir = directory_path

    print(results_dir)

    # Data frame to store the summarized results for saturation over all categories
    summarized_data, summarized_data_correction = summarize_saturated_results_categories(directory_path, dataset_name)
    
    if not summarized_data.empty:  
        # Rearrange summarized data
        rearranged_data = rearrange_dataframe(summarized_data)

        # Save rearranged data
        csv_file_path = os.path.join(results_dir, f"{dataset_name}_summary_saturation.csv")
        # Write DataFrame to CSV file
        rearranged_data.to_csv(csv_file_path, index=False)

        #Map rearranged data to newick string 
        if newick_string:
            # map rearranged data to it
            #print(newick_string)
            result_nexus_file = os.path.join(results_dir, f"{dataset_name}_summary_saturation_tree.nex")
            write_nexus_file(newick_string,result_nexus_file,rearranged_data)
    else:
        print("No saturated branches for ", dataset_name, "!")

    if not summarized_data_correction.empty:  
        # Rearrange summarized data
        rearranged_data = rearrange_dataframe_correction(summarized_data_correction)

        # Save rearranged data
        csv_file_path = os.path.join(results_dir, f"{dataset_name}_summary_saturation_bonferroni_corrected.csv")
        # Write DataFrame to CSV file
        rearranged_data.to_csv(csv_file_path, index=False)

        #Map rearranged data to newick string 
        if newick_string:
            # map rearranged data to it
            #print(newick_string)
            result_nexus_file = os.path.join(results_dir, f"{dataset_name}_summary_saturation_bonferroni_corrrected.nex")
            write_nexus_file(newick_string,result_nexus_file,rearranged_data)
        
    return summarized_data

def process_gene_directory_z_scores(directory_path, dataset_name, newick_string=None):

    # Data frame to store tz scores over all categories
    summarized_z_score_data = summarize_z_scores_categories(directory_path, dataset_name)
    
    if not summarized_z_score_data.empty:  
        # Save data
        csv_file_path = os.path.join(directory_path, f"{dataset_name}_summarized_z_scores.csv")
        # Write DataFrame to CSV file
        summarized_z_score_data.to_csv(csv_file_path, index=False)

    #     #Map rearranged data to newick string 
    #     if newick_string:
    #         # map rearranged data to it
    #         #print(newick_string)
    #         result_nexus_file = os.path.join(directory_path, f"{dataset_name}_summary_z_tree.nex")
    #         write_nexus_file(newick_string,result_nexus_file,rearranged_data)
        
    return summarized_z_score_data

def process_gene_directory_components(directory_path, dataset_name, newick_string=None):

    # Data frame to store data over all categories
    summarized_data = summarize_components_categories(directory_path, dataset_name)
    
    # if not summarized_data.empty:  
    #     # Save data
    #     csv_file_path = os.path.join(directory_path, f"{dataset_name}_summarized_components.csv")
    #     # Write DataFrame to CSV file
    #     summarized_data.to_csv(csv_file_path, index=False)

    #     #Map rearranged data to newick string 
    #     if newick_string:
    #         # map rearranged data to it
    #         #print(newick_string)
    #         result_nexus_file = os.path.join(directory_path, f"{dataset_name}_summary_z_tree.nex")
    #         write_nexus_file(newick_string,result_nexus_file,rearranged_data)
        
    return summarized_data

if __name__ == "__main__":

    """ Set directories and input files """
    # Get the current working directory
    current_directory = os.getcwd()

    # Get directory of the considered example: SSU rRNA tree
    ssu_msa_dir = os.path.join(current_directory, "../SSU_MSA/")

    # Specify the path to results directories
    results_dir = os.path.join(ssu_msa_dir, "results")
    # Name the dataset
    dataset_name = "SSU-MSA-Suppl3-renamed.fasta"

    print(summarize_results_categories(results_dir,dataset_name))