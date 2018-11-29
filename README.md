# A python script for branch-site positive-selection analysis

## Requirement
* Python >= 2.7 (3.x is recommended)
* ETE3 package

## Usage
You can get a summary of available command-line options with `branch_site_pipeline.py -h`
```
$ python branch_site_pipeline.py -h
usage: branch_site_pipeline.py [-h] [-i INPUT_FILE] [-o OUTPUT_DIR] [-v]
                               [-l LOGFILE] [-t THREADS]
[...]
```
### Input
You need to prepare a text file including these information:
1. cDNA alignment file (*.fasta, *.paml);
2. tree file (newick format);
3. Any specific nodes in the tree. (Cautions: all nodes need to be monophyletic). If you want to test all branches, just leave blank here.

Note: `test.txt` in the source package is an example input file.

### Output
Two files will be output in the output directory:
1. `all_results.txt` including all calculation results and nodes diagrams;
2. `positive_selection_tree.png` showing all nodes with positive selection.