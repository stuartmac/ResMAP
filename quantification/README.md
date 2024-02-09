# ResMAP Library Quantification

The python script `count-sequences.py` is used to count the occurence of each mutant in the library. The script is run as follows:

```bash
# First create a virtual environment and install the required packages
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Run the script
python count-sequences.py -i <input_file> -o <output_file> --library <library_file>
```

The input file is a SAM or BAM file of the mapped sequencing library. The library file should be a csv file containing the library sequences together with an identifier.

A sample of the PfKRS sequencing data is provided in the `data` directory.

To run the script on the sample data, use the following command:

```bash
python count-sequences.py -i data/PkmutWTd1_1pc.sam -o data/PkmutWTd1_1pc-counts.tsv --library data/pfkrs-n20-library.csv
```

The output file will be `PkmutWTd1_1pc-counts.tsv` and will contain the counts of each mutant in the library.