# ResMAP

This repository contains the code for the analysis of the deep mutational screening data from the ResMAP project as detailed in the publication:

 - Richard J. Wall, Stuart A. MacGowan, Irene Hallyburton, Sowmya Ajay Castro, Gourav Dey, Rachel Milne, Stephen Patterson, Jody Phelan, Aisha Syed, Natalie Wiedemar and Susan Wyllie. *ResMAP – a saturation mutagenesis platform enabling parallel profiling of target-specific resistance conferring mutations in Plasmodium.* Submitted.

ResMAP aims to identify mutations in Plasmodium that confer resistance to antimalarial drugs, providing insights critical for developing novel therapeutics.

## Repository Structure

The `quantification`, `analysis`, and `data` directories in this repository encompass the entire workflow from data quantification to analysis:
- `quantification`: Contains `count-sequences.py` for library quantification.
- `analysis`: R Markdown notebooks for data analysis.
- `data`: Includes raw and processed data files.

## Data

The `data` directory structure is as follows:
- `PfKRS-mutlib-final-calculations-*.tsv`: Library counts for various conditions.
- `PkmutWTd1_1pc.sam`: Sample sequencing data.
- `pfkrs-n20-library.csv`: Library design file.

## Quantification

`count-sequences.py` in the `quantification` directory is used for counting mutants. See the [quantification README](quantification/README.md) for usage details.

## Analysis

R Markdown notebooks in the `analysis` directory outline the analysis workflow. Numbered sequentially, these notebooks are detailed in the [analysis README](analysis/README.md), with rendered versions available for viewing.

## License

This repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
```
