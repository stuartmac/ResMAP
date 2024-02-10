# ResMAP

This repository contains the code for the analysis of the deep mutational screening data from the ResMAP project reported in:

"ResMAP – a saturation mutagenesis platform enabling parallel profiling of target-specific resistance conferring mutations in Plasmodium"

by Richard J. Wall, Stuart A. MacGowan, Irene Hallyburton, Sowmya Ajay Castro, Gourav Dey, Rachel Milne, Stephen Patterson, Jody Phelan, Aisha Syed, Natalie Wiedemar and Susan Wyllie

ResMAP is a platform designed to systematically identify mutations in Plasmodium that confer resistance to antimalarial drugs, offering critical insights for the development of novel therapeutics.

## Repository Structure

The repository is organised into three main directories: `quantification`, `analysis` and `data`. The `quantification` directory contains the code for library quantification, the `analysis` directory contains the code for the analysis of the deep mutational screening data, and the `data` directory contains the library summary data from the deep mutational screening.

## Data

data/
├── PfKRS-mutlib-final-calculations-706_231117.tsv  # library counts for compound 706
├── PfKRS-mutlib-final-calculations-clado_231117.tsv  # library counts for cladosporin
├── PfKRS-mutlib-final-calculations-wildtype_231117.tsv  # library counts pre-selection
├── PkmutWTd1_1pc-counts.tsv
├── PkmutWTd1_1pc.sam  # sample sequencing data
├── PkmutWTd1_1pc.sam.bai
├── codon-table.tsv  # resource for Rmd analysis
├── mutlib-1187736dfa7e94e0c60d8ef17aa56979.tsv
└── pfkrs-n20-library.csv  # library design file for quantification

## Quantification

The python script for library quantification is located in the `quantification` directory. The script `count-sequences.py` is used to count the occurence of each mutant in the library. Further details on the script and its usage can be found in the [README](quantification/README.md) in the `quantification` directory.

## Analysis 

The R Markdown notebooks for the analysis of the deep mutational screening data are located in the `analysis` directory. The notebooks are numbered in the order they should be run, and further details on each notebook can be found in the [README](analysis/README.md) in the `analysis` directory.

Rendered versions of the notebooks can also be found in the `analysis` directory. These notebooks contain the results of the analysis and can be viewed directly in the browser through GitHub.

# License

This project code is open source and available under the [MIT License](LICENSE).
[]: # (end)
```