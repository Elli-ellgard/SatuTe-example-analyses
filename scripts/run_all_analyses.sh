#!/bin/bash

for dir in ../example/results*; do
  if [ -d "$dir" ]; then
    rm -r "$dir"
  fi
done

python=python3

# branch specific sliding analysis
cd branch_specific_sliding_window_analysis/
$python  script_sliding_window_analysis.py
$python  script_visualise_results.py 

# per category analyses
cd ../per_category_analysis/
$python script_category_analysis.py 

# per region analysis
cd ../per_alignment_region_analysis
$python script_per_region_analysis.py

#missing separate analysis

# z-score differences between  branches
cd  ../zscore_differences_between_branches
$python script_per_region_zscore_differences_between_branches.py

# z-score differnces between topologies
cd  ../zscore_differences_between_topologies
$python script_per_region_zscore_differences_between_topologies.py