# Utility-Oriented-String-Indexing
This repository contains the source code of the paper titled "Utility-Oriented String Indexing". The main problem in the paper is the USI problem, which focuses on preprocessing a string with numerical utilities into a compact data structure, enabling efficient computation of the global utility of query patterns based on their occurrences and associated utilities. The data structure we propose requires also finding top-K frequent substrings from a long string.

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

### Usage for ET and AT:
```bash
./et <text_file> <K>
./at <text_file> <K> <s>
```

### Example commands for ET and AT:
```bash
./et ../data/string ../data/weight_string 3
./at ../data/string ../data/weight_string 3 2
```


## 2. USI Problem
```bash
cd ../usi
```
### Input parameters:
- `<text_file>` : Input string file.
- `<weight_file>` : Weight file for the string.
- `<K>` : The number of top frequent substrings to find.
- `<s>` : The number of partitions.
- `<pattern file>` : Query pattern file.

  
### Usage for UET, UAT, BSL1, BSL2, BSL3, and BSL4:
```bash
./uet <text_file> <weight_file> <K> <pattern file>
./uat <text_file> <weight_file> <K> <s> <pattern file>
./bsl1 <text_file> <weight_file> <pattern file>
./bsl2 <text_file> <weight_file> <K> <pattern file>
./bsl3 <text_file> <weight_file> <K> <pattern file>
./bsl4 <text_file> <weight_file> <K> <pattern file>
```

### Example Commands UET, UAT, BSL1, BSL2, BSL3, and BSL4:

```bash
./uet ../data/string ../data/weight_string 3 ../data/pattern
./uat ../data/string ../data/weight_string 3 2 ../data/pattern
./bsl1 ../data/string ../data/weight_string ../data/pattern
./bsl2 ../data/string ../data/weight_string 3 ../data/pattern
./bsl3 ../data/string ../data/weight_string 3 ../data/pattern
./bsl4 ../data/string ../data/weight_string 3 ../data/pattern
```


