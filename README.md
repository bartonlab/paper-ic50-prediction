# **Overview**  

This repository contains the necessary data files, source code, and notebooks required to generate the figures included in the manuscript, *Predicting viral sensitivity to antibodies using genetic sequences and antibody similarities* ([PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014095)).  

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

This study is published from [PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014095).  

---

## **Repository Contents**

This repository contains the following directories:

- **`raw-data/`** – Raw input data retrieved from LANL databases
- **`processed-data/`** – Processed data files used for simulation
- **`src/`** – Julia source code files
- **`note/`** – Jupyter notebooks for data processing, simulation, and visualization

### **`raw-data/`**
This directory contains the raw data used in this project, retrieved from LANL databases. These include:

- neutralization data from the **CATNAP database**
- HIV-1 envelope sequence data from the **LANL HIV database**

These raw datasets are used as the starting point for downstream preprocessing and simulation.

### **`processed-data/`**
This directory contains the processed data files used for simulation.

These files are generated from the raw data using the following notebooks:

- **`process_non-alignment-features.ipynb`** – processes the non-alignment-based features
- **`process.ipynb`** – processes the remaining required data files

### **`GNL.ipynb`**
This notebook runs the neutralization matrix imputation simulations.

It imputes missing neutralization values and also evaluates imputation performance by randomly withholding observed values, imputing them, and exporting the imputed results for comparison with the true withheld values.

This notebook also analyzes the rank dependency of imputation performance.

### **`simple_visualization_of_imputation_results.ipynb`**
This notebook visualizes the imputation results.

It includes:

- scatter plots comparing imputed values with the true withheld values
- figures showing the rank dependency of imputation accuracy

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

