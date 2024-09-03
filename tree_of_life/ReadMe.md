# Tree of Life Analyses

Ensure that all prerequisites are installed before proceeding. You can either run the Tree of Life analyses individually or you can execute all analyses at once. To execute all scripts at once, navigate to `scripts` and use the following command:

```bash
./run-all-ToL-analyses.sh
```

If you prefer to run the analyses individually, navigate to the specified directory (location) and execute the script. The resulting figures will closely resemble those presented in the manuscript.

## 1. Branch-Specific Per-Category Analysis

Perform per-category analysis for the Eukaryota and Yeast branches in the 2D or 3D Tree of Life as follows:

### 1.1. 16S rRNA-Based 3D ToL

1. **Script:** `pca-rna-based-3D-ToL.py`
2. **Location:** `scripts/branch_specific_per_category_analysis`
3. **Run:**

    ```bash
    python3 pca-rna-based-3D-ToL.py
    ```

4. **Output Directory:** `results_per_category_analysis/rRNA_based_3D_tree`

### 1.2. Protein-Based 2D ToL

1. **Script:** `pca-protein-based-2D-ToL.py`
2. **Location:** `scripts/branch_specific_per_category_analysis`
3. **Run:**

    ```bash
    python3 pca-protein-based-2D-ToL.py
    ```

4. **Output Directory:** `results_per_category_analysis/protein_based_2D_tree`

## 2. Branch-Specific Sliding-Window Analysis

Perform sliding-window analysis for the Eukaryota and Yeast branches in the 2D or 3D Tree of Life as follows:

### 2.1. 16S rRNA-Based 3D ToL

1. **Script:** `swa-rna-based-3D-ToL.py`
2. **Location:** `scripts/branch_specific_sliding_window_analysis`
3. **Adjust Parameters:** Modify window size and branch selection as needed.
4. **Run:**

    ```bash
    python3 swa-rna-based-3D-ToL.py
    ```

5. **Output Directory:** `results_sliding_window_analysis/rRNA_based_3D_tree`

### 2.2. Protein-Based 2D ToL

1. **Script:** `swa-protein-based-2D-ToL.py`
2. **Location:** `scripts/branch_specific_sliding_window_analysis`
3. **Adjust Parameters:** Modify window size and branch selection as needed.
4. **Run:**

    ```bash
    python3 pca-protein-based-2D-ToL.py
    ```

5. **Output Directory:** `results_sliding_window_analysis/protein_based_2D_tree`

## 3. Z-Score Differences between Branches

Calculate z-score differences between the Eukaryota and Yeast branch in the 2D Tree of Life as follows:

### 3.1. Protein-Based 2D ToL

1. **Script:** `zscore-differences-Eukaryota-Yeast-in-protein-based-2D-ToL.py`
2. **Location:** `scripts/zscore_differences_between_branches`
3. **Run:**

    ```bash
    python3 zscore-differences-Eukaryota-Yeast-in-protein-based-2D-ToL.py
    ```

4. **Output Directory:** `results_zscore_differences_branches/protein_based_2D_tree`

## 4. Z-Score Differences between Topologies

Calculate the z-score differences for the Eukaryota branch in a 2D or 3D tree as follows:

### 4.1. 16S rRNA-Based 3D ToL

1. **Script:** `zscore-differences-rna-based-2D-vs-3D-ToL.py`
2. **Location:** `scripts/zscore_differences_between_toplogies/`
3. **Run:**

    ```bash
    python3 zscore-differences-rna-based-2D-vs-3D-ToL.py
    ```

4. **Output Directory:** `results_zscore_differences_topologies/rRNA_based_3D_tree`

### 4.2. Protein-Based 2D ToL

1. **Script:** `zscore-differences-protein-based-3D-vs-2D-ToL.py`
2. **Location:** `scripts/zscore_differences_between_toplogies/`
3. **Run:**

    ```bash
    python3 zscore-differences-protein-based-3D-vs-2D-ToL.py
    ```

4. **Output Directory:** `results_zscore_differences_topologies/protein_based_2D_tree`