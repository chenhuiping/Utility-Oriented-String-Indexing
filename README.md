# Utility-Oriented String Indexing
This repository contains the source code of the paper titled "Indexing Strings with Utilities". The main problem in the paper is the USI problem, which focuses on preprocessing a string with numerical utilities into a compact data structure, enabling efficient computation of the global utility of query patterns based on their occurrences and associated utilities. The data structure we propose requires also finding top-K frequent substrings from a long string.

There is source code for:

1. Top-K Frequent Substring Mining: Algorithms ET, AT, TT, and SH.

2. USI Problem: Data structures UET, UAT, BSL1, BSL2, BSL3, and BSL4.

# Requirements
- A GNU/Linux system
- A modern C++17 ready compiler
- This project requires the installation of the [SDSL-lite](https://github.com/simongog/sdsl-lite) library. Please follow these steps to set it up:
```bash
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
```

After installing SDSL-lite, Update the configuration file located at `path.cfg` to set the correct library path.

For example: `SDSL_LIB_DIR=/path/to/sdsl-lite` 


# Installation
```bash
./compile.sh
```


# How to run

## 1. Top-K Frequent Substring Mining
```bash
cd topk
```
### Input parameters:
- `<text_file>` : Input string file.
- `<K>` : The number of top frequent substrings to find.
- `<s>` : The number of partitions.

### Usage for ET, AT and SH:
```bash
./et <text_file> <K>
./at <text_file> <K> <s>
./sh <text_file> <K>
```

### Example commands for ET, AT and SH:
```bash
./et ../data/example_string 3
./at ../data/example_string 3 2
./sh ../data/example_string 3 
```
### TT algorithm:
The TT algorithm, as implemented in the top-k-compress project focuses on identifying the top-k frequent patterns in data streams. For a comprehensive understanding of the TT algorithm's methodology and its applications, please refer to the [topk-k-compress](https://github.com/pdinklag/top-k-compress) GitHub Repository. We utilised the `topk-lz78` algorithm to conduct our experiment.

## 2. USI Problem
```bash
cd ../usi
```
### Input parameters:
- `<text_file>` : Input string file.
- `<weight_file>` : Weight file for the string.
- `<K>` : The number of top frequent substrings to find.
- `<s>` : The number of partitions.
- `<pattern_file>` : Query pattern file.

  
### Usage for UET, UAT, BSL1, BSL2, BSL3, and BSL4:
```bash
./uet <text_file> <weight_file> <K> <pattern_file>
./uat <text_file> <weight_file> <K> <s> <pattern_file>
./bsl1 <text_file> <weight_file> <pattern_file>
./bsl2 <text_file> <weight_file> <K> <pattern_file>
./bsl3 <text_file> <weight_file> <K> <pattern_file>
./bsl4 <text_file> <weight_file> <K> <pattern_file>
```

### Example Commands UET, UAT, BSL1, BSL2, BSL3, and BSL4:

```bash
./uet ../data/example_string ../data/example_weight 3 ../data/example_pattern
./uat ../data/example_string ../data/example_weight 3 2 ../data/example_pattern
./bsl1 ../data/example_string ../data/example_weight ../data/example_pattern
./bsl2 ../data/example_string ../data/example_weight 3 ../data/example_pattern
./bsl3 ../data/example_string ../data/example_weight 3 ../data/example_pattern
./bsl4 ../data/example_string ../data/example_weight 3 ../data/example_pattern
```
# Contact

For questions or support, please contact:  
Huiping Chen (<h.chen.13@bham.ac.uk>), Roberto Grossi (<roberto.grossi@unipi.it>), Veronica Guerrini (<veronica.guerrini@unipi.it>), and Solon P. Pissis (<solon.pissis@cwi.nl>).


