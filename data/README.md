# Data Description

This directory contains the datasets used in this project. Below is an explanation of the structure of each input file and the data sources.

# Supported File Formats

**1. Input String File `<text_file>`:** contains a single long string of characters.

Example:
```bash
banana
```
**2. Weight File `<weight_file>`:** provides a numerical weight for each character in the corresponding string file, weights are separated by spaces and must align with the string length.

Example:
```bash
  1 3 4 3 4 1
```
**3. Query Pattern File `<pattern file>`:** lists query patterns (substrings) to be used for utility-based retrieval, each line represents a substring.

Example:
```bash
ba
a
ana
n
an
```

## Real Datasets

The following real-world datasets are included or referenced in this project:

- [**ADV**](https://raw.githubusercontent.com/UmerTariq1/CTR-and-Clicks-Prediction-Using-Text-and-Numerical-Data/master/Dataset.csv)  : Included in the folder as `ADV.txt`: Dataset for predicting click-through rates (CTR) using text and numerical data.

- [**IOT**](https://ieee-dataport.org/open-access/crawdad-unmblebeacon): Captures signals from Internet of Things (IoT) devices.

- [**XML**](https://pizzachili.dcc.uchile.cl/texts.html): XML data used for text compression and indexing benchmarks.

- [**HUM**](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/): Genomic data for human sequences.

- [**ECOLI**](https://www.ebi.ac.uk/ena/browser/view/ERR022075): Genomic data for *Escherichia coli* sequences.


