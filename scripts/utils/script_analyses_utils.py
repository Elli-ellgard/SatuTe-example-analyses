import os
import subprocess


   
def find_fasta_file(gene_name, source_folder):
    # Iterate over files in the source folder
    for filename in os.listdir(source_folder):
        if gene_name in filename and filename.endswith(".fasta"):
            return os.path.join(source_folder, filename)
    return None

def run_external_command(command_args):
    subprocess.run(command_args)

def find_tree_file(gene_name, source_folder):
    # Iterate over files in the source folder
    for filename in os.listdir(source_folder):
        if gene_name in filename and filename.endswith(".treefile"):
            return os.path.join(source_folder, filename)
    return None

def find_file_with_suffix(gene_name, suffix, source_folder):
    # Iterate over files in the source folder
    for filename in os.listdir(source_folder):
        if gene_name in filename and filename.endswith(suffix):
            return os.path.join(source_folder, filename)
    return None


def find_file_with_suffix_in_directory(suffix, source_folder):
    # Iterate over files in the source folder
    for filename in os.listdir(source_folder):
        if filename.endswith(suffix):
            return os.path.join(source_folder, filename)
    return None


def get_regions_from_annotation(annotation_file):
    # Read the annotation file to get site and region_name information
    region_info = {}
    with open(annotation_file) as f:
        next(f)  # Skip header
        for line in f:
            site, region_name = line.strip().split(",")
            site = int(site)
            region_name = region_name.strip('"')
            region_info[region_name] = region_info.get(region_name, []) + [site]
    return region_info