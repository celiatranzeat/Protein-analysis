# Protein Analysis Project

**Author:** Célia TRANZEAT  
**Date:** 24 October 2025  
**Contact:** celia.tranzeat-lecourtois@etu.univ-cotedazur.fr  

---

## What you need to use the application ?

This application requires a working **Python 3** environment with the following modules installed:

- `biopython`  
- `ipwidgets`
- `matplotlib` 
-  `Projet_Celia_TRANZEAT.py` 
- A Jupyter Notebook environment (such as Anaconda)

---

## How to make it run ?

To run the application, follow these steps:

1. Launch **Jupyter Notebook**.  
2. Open the notebook file **`Interface_graphique_Celia_TRANZEAT.ipynb`**.  
3. In the interface graphique cell, run the code by pressing **Shift + Enter**.  
4. When prompted, **enter a valid protein ID** (for example, a UniProt or NCBI ID such as `O23729`).  
5. The interface will then display various options and graphical outputs automatically.

---

## How it works ?

The project is based on a Python module named **`Projet_Celia_TRANZEAT.py`**, which contains all the functions used by the graphical interface.  

### Input  
- A **protein identifier** .  
- The program fetches the corresponding protein record automatically using **ExPASy** and **BioPython**.  

### Process  
Once the record is retrieved, the program performs several analyses:
- Extraction of key biological information:  
  *Uniprot ID, gene name, accessions, taxonomy, sequence length, comments, and subcellular location.*
- Calculation of biochemical properties such as:  
  - **Amino acid composition** (percentage of each residue)  
  - **Hydrophobicity profile** (Kyte & Doolittle scale)  
  - **Net charge as a function of pH** (0–14 range)  
  - **Detection of N-glycosylation motifs** (N-X-[S/T], with X ≠ P)

### Output  
The program generates:  
- A **dictionary** containing all extracted information.  
- A **bar chart** showing amino acid composition.  
- A **hydrophobicity curve** along the sequence.  
- A **net charge vs pH graph**.  
- A **list of N-glycosylation motif positions** detected in the protein.  

All plots are displayed within the Jupyter interface using **matplotlib**.  

---

## Licence

This project is distributed for educational and academic use under a  
**Creative Commons Attribution-NonCommercial (CC BY-NC)** license.  

You are free to use and modify the code for non-commercial purposes,  
provided that proper credit is given to the author.  

---

## Citation

If you use this application or any part of its code in your work, please cite the following reference:

> Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009).  
> **Biopython: freely available Python tools for computational molecular biology and bioinformatics.**  
> *Bioinformatics*, 25(11), 1422–1423.  
> doi:[10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)
