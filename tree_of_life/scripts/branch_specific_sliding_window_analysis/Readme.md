# Tree of Life Analyses

## 1. Branch-Specific Sliding Window Analysis

For performing branch-specific sliding window analysis on 2D or 3D Tree of Life (ToL) data, follow the steps below:

### 1.1. 16S rRNA-Based 3D ToL

**Description:**


**How to Run:**

1. **File to Run:** `swa-rna-based-3D-ToL.py`
2. Adjust the script parameters as needed, including window size and branch selection.
3. Execute the script by running the following command in your terminal or command prompt:

    ```bash
    python3 swa-rna-based-3D-ToL.py
    ```

**Output:**

- All output files will be stored in the directory `../results_sliding_window_analysis/rRNA_based_3D_tree`.
- To generate the original figure from the paper, execute:

    ```bash
    Rscript <your_script_name>.R
    ```

### 1.2. Protein-Based 2D ToL

**Description:**


**How to Run:**

1. **File to Run:** `swa-protein-based-2D-ToL.py`
2. Adjust the script parameters as needed, including window size and branch selection.
3. Execute the script by running the following command in your terminal or command prompt:

    ```bash
    python3 swa-protein-based-2D-ToL.py
    ```

**Output:**

- All output files will be stored in the directory `../results_sliding_window_analysis/protein_based_2D_tree`.
- To generate the original figure from the paper, execute:

    ```bash
    Rscript <your_script_name>.R
    ```



