# **Overview**  

This repository contains the necessary data files, source code, and notebooks required to generate the figures included in the manuscript, *Predicting viral sensitivity to antibodies using genetic sequences and antibody similarities* ([bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.08.669352v1.full)).  

Our proposed method, **Grouped Neutralization Learning (GNL)**, predicts neutralization values from partial neutralization data (e.g., IC50, IC80, etc.) and viral surface protein sequences. In this study, we focus specifically on HIV-1 sequences and anti-HIV-1 titer values.  

<!-- The primary purpose of this repository is to reproduce the results reported in the manuscript. Therefore, we provide notebooks and source code that perform analysis and evaluate the method’s performance. However, we have a separate repository dedicated solely to predicting neutralization values given new neutralization data and viral sequences. If you are primarily interested in the application of our software, please visit the repository for *Title* ([bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.08.669352v1.full)).  i-->

---

## **Manuscript Information**  

**Predicting viral sensitivity to antibodies using genetic sequences and antibody similarities**  
Kai S. Shimagaki<sup>1,2</sup>, Gargi Kher<sup>1</sup>, Rebecca M. Lynch<sup>3</sup>, and John P. Barton<sup>1,2,#</sup>  

<sup>1</sup> Department of Computational and Systems Biology, University of Pittsburgh School of Medicine, USA.  
<sup>2</sup> Department of Physics and Astronomy, University of Pittsburgh, USA.
<sup>3</sup> Department of Microbiology, Immunology and Tropical Medicine, School of Medicine and Health Sciences, George Washington University, USA.
<sup>#</sup> Correspondence: [jpbarton@pitt.edu](mailto:jpbarton@pitt.edu)  

The preprint is available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.08.669352v1.full).  

---

## **Repository Contents**  

This repository contains the following directories:  

- **`data/`** – Contains all data files  
  - **`data/raw/`** – Raw data before processing  
  - **`data/processed/`** – Processed data for analysis  
- **`src/`** – Source code files  
- **`notebooks/`** – Jupyter notebooks for data processing, analysis, and figure generation  
- **`figures/`** – Generated figures, which can be produced using `figures/figures.ipynb`  

The **`data/raw/`** directory contains neutralization data obtained from the **CATNAP database** (retrieved in November 2024). It also includes HIV-1 surface protein sequences obtained from the **HIV database** and aligned using the **HIV-Align tool** provided by LANL.  

---

## **Notebooks Overview**  

The **`notebooks/`** directory contains three primary Jupyter notebooks:  

### 1. **`process_data.ipynb`**  
- Preprocesses raw neutralization data from the **CATNAP database** (LANL) and viral sequences in **FASTA format**.  
- Generates multiple processed files necessary for neutralization imputation and subsequent analysis.  
- Outputs include:  
  - **One-hot encoded sequences** (filtered for data quality)  
  - **Projected sequences** (dimensionality reduction of one-hot encoded sequences)  
  - **PCA space representation** (projected sequences onto PCA space)  
  - **Neutralization matrices** (rows: antibodies, columns: viruses) with and without filtering  
  - **Accessory files** for further analysis  

### 2. **`analysis.ipynb`**  
- Evaluates imputation performance by withholding neutralization values **at random** at multiple withholding rates (effectively 1 - observed data fraction).  
- Partitions imputed data based on:  
  - Number of observations per entry  
  - Standard deviation (quantifying heterogeneity of neutralization values, a key parameter for assessing imputation performance)  

### 3. **`figures.ipynb`**  
- Uses results from `analysis.ipynb` to generate **visualizations** of the imputation analysis.  
- Outputs figures stored in the **`figures/`** directory.  

---

## **Software Dependencies**  

The code is implemented in **Julia (version 1.8)** and uses the following standard libraries:  
- **Random**, **LinearAlgebra**, **StatsBase**, **Distributed**, **Statistics**, **DelimitedFiles**, **Printf**

---

## **License**  

This repository is **dual-licensed** under:  
- **[GPL-3.0](LICENSE-GPL)** for source code  
- **[CC0 1.0](LICENSE-CC0)** for figures, documentation, and presentation of data  

---

### **Additional Notes**  
- If you encounter any issues or have questions, please reach out via **GitHub issues** or email the corresponding author.  
<!-- - For users interested in applying the method to their own data, please refer to the separate repository: [**Predicting viral sensitivity to antibodies using genetic sequences and antibody similarities** (Link)].  -->

