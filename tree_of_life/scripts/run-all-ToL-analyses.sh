#!/bin/bash

for dir in ../results*; do
  if [ -d "$dir" ]; then
    rm -r "$dir"
  fi
done

python=python3

# branch specific per category analyses
cd branch_specific_per_category_analysis/
$python  pca-protein-based-2D-ToL.py
$python  pca-rna-based-3D-ToL.py

# branch specific sliding analysis
cd ../branch_specific_sliding_window_analysis/
$python  swa-protein-based-2D-ToL.py
$python  swa-rna-based-3D-ToL.py

# z-score differences between  branches
cd  ../zscore_differences_between_branches
$python zscore-differences-Eukaryota-Yeast-in-protein-based-2D-ToL.py

# z-score differnces between topologies
cd  ../zscore_differences_between_topologies
$python zscore-differences-rna-based-2D-vs-3D-ToL.py
$python zscore-differences-protein-based-3D-vs-2D-ToL.py