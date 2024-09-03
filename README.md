# SatuTe Example Analyses

[SatuTe](https://github.com/Elli-ellgard/SatuTe) (Saturation Test) is a Python-based tool designed to evaluate the presence of phylogenetic information in analyses. Saturation occurs when multiple substitutions obscure true genetic distances, potentially leading to artifacts and errors in phylogenetic inference. SatuTe introduces a new measure that extends the concept of saturation between two sequences to a theory of saturation between subtrees. The test implemented in SatuTe quantifies whether a given alignment provides sufficient phylogenetic information shared between two subtrees connected by a branch in a phylogeny.

Using the output from SatuTe, you can perform various downstream analyses to gain deeper insights into the phylogenetic signal by addressing different questions.

## 1. Repository Structure

The repository is organized as follows:

- **/example/**: Contains small example datasets and the output generated by SatuTe. This folder allows you to follow along with the different types of analyses.

- **/scripts/**: Includes all scripts required for performing the various analyses, such as per-category, sliding-window, and per-alignment-region analyses. Each type of analysis has its own subfolder, and the scripts can be run to generate examples within the `/example/` folder. Note that installation of SatuTe and IQ-TREE is not necessary to run these examples, as the required outputs are already provided.

- **/tree_of_life/**: Contains the data and scripts used to generate the outputs presented in the associated paper. This includes detailed instructions and resources for replicating the findings.

## 2. Prerequisites

Before running any scripts, ensure you have the following software installed:

- **Python 3.x**: Confirm that Python is installed on your system.

### Additional Tools for Tree of Life Analysis

If you are planning to run the Tree of Life analysis, you'll need these additional tools:

- **SatuTe**
- **IQ-TREE**

#### Install SatuTe using pipx

1. **Install pipx**: If you don't have pipx installed, you can install it using pip:

   ```bash
   pip install pipx
   ```

2. **Ensure pipx is set up correctly**:

   ```bash
   pipx ensurepath
   ```

3. **Install SatuTe using pipx**: Once pipx is installed, you can use it to install SatuTe:

   ```bash
   pipx install satute
   ```

4. **Test the installation**: After installation, verify that SatuTe is installed correctly by checking its version:

   ```bash
   satute --version
   ```

For more detailed instructions and information about pipx, refer to the [official pipx documentation](https://pypa.github.io/pipx/).

#### IQ-TREE Installation

1. **Download IQ-TREE** from the [official website](http://www.iqtree.org/#download).
  
2. **Follow the installation instructions** provided on the website for your operating system.

3. **Test the IQ-TREE installation**: After installing IQ-TREE, verify the installation by checking its version:

   ```bash
   iqtree2 --version
   ```

## 3. Types of Analyses

Using the output from SatuTe, you can perform various downstream analyses:

### 3.1. Per-Category Analysis

If an evolutionary model with rate heterogeneity is used, each site is assigned to the rate category with the highest posterior probability. For each category $c$, SatuTe applies the test for phylogenetic information on the rescaled phylogenetic tree and the subalignment of the considered category. During this process, SatuTe calculates the variance estimator $\hat{\sigma}^2_{1,c}$ for each category $c$, enabling you to compare the phylogenetic signals present in each rate category.

### 3.2. Branch-Specific Sliding-Window Analysis

SatuTe also supports branch-specific analyses. To gain a more detailed understanding of changes in phylogenetic information, you can perform a sliding-window analysis with a specific window size. This approach is effective in detecting a minority of sites affected by saturation.

### 3.3. Per-Alignment Region Analysis

When the alignment is composite —such as a concatenation of different genes, proteins, or other partitions— a key question is whether the selected alignment regions are phylogenetically informative within the reconstructed tree topology. A per-alignment-region analysis can help address this question.

### 3.4. Z-Score Differences

By comparing the z-scores obtained from different branches, you can identify potential information loss and examine the differences. For instance, you might explore per-region z-score differences between an external branch and an internal branch. Beyond branch comparison, z-score differences can also help determine whether each region supports one of two given topologies.
