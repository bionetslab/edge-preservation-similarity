# Edge-Preservation-Similarity
This git provides an exact and an approximated algorithm for computing the edge-preservation similarity between rooted, unordered, node-labeled trees.

## GUROBI










# DIGEST
The source code for [DIGEST](https://digest-validation.net/) (validation of **di**sease and **ge**ne **s**ets or clus**t**erings). It greatly facilitates in silico validation of gene and disease sets or clusterings via fully automated validation pipelines comprising disease and gene ID mapping, enrichment
analysis, comparisons of shared genes and variants, and background distribution estimation. Moreover, functionality is provided to automatically update the external databases used by the pipelines.

## Setup
1. Install git
```
pip install git
```
2. Clone this repository
```
git clone git@github.com:digest-env/digest.git
```
3. Setup enviroment

3.1. Import yml file to enviroment
```
conda env create -f environment.yml
```
3.2. Setup manually

3.2.1 Setup enviroment
```
conda create --name digest python==3.8
conda activate digest
```
3.2.2. Install dependancies
```
pip install pandas numpy scipy seaborn biothings_client gseapy
python -m pip install psutil
```
## Setup files
To make sure, that all mappings are up to date, run the setup script. This will retrieve the mappings from the api. Runtime: ~1 Minute. This is recommended as the files on the api will be kept updated.
```
python3 setup.py
```
Alternatively you can setup the files while creating them from scratch. This is not recommended, as depending on the system, it could run up to 3 hour.
```
python3 setup.py -s="create"
```
## Run DIGEST
### Run in terminal
```
usage: python3 single_validation.py [required arguments] [optional arguments]
required arguments:
  -r REFERENCE, --reference REFERENCE
                        [Only for mode set-set] Reference file or id. 
  -ri REFERENCE_ID_TYPE, --reference_id_type REFERENCE_ID_TYPE
                        [Only for mode set-set] Reference id type. See possible options below.
  -t TARGET, --target TARGET
                        Target file with set or clusters.
  -ti TARGET_ID_TYPE, --target_id_type TARGET_ID_TYPE
                        Target id type. See possible options below.
  -m {set,set-set,cluster}, --mode {set,set-set,cluster}
                        Desired mode. See possible options below.
optional arguments:
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -dg {jaccard,overlap}, --distance_measure {jaccard,overlap}
                        Distance measure. [Default=jaccard]
  -e, --enriched        Set flag, if only enriched attributes of the reference should be used.
  -c RUNS, --runs RUNS  Number of runs with random target values for p-value calculation.
  -b {complete,term-pres}, --background_model {complete,term-pres}
                        Model defining how random values should be picked. See possible options below.
  -pr REPLACE, --replace REPLACE
                        Percentage of how many of the original ids should be replaced with random ids. [Default=100]
  -v, --verbose         Set flag, if additional info like ids without assigned attributes should be printed.
  -p, --plot            Set flag, if plots should be created.
  -h, --help            show this help message and exit
----------------------------------------------------------------------------
supported id types
  for genes		entrez, ensembl, symbol, uniprot
  for diseases		mondo, omim, snomedct, umls, orpha, mesh, doid, ICD-10
supported modes
  set			Compare similarity inside the set. Either genes or diseases.
  set-set		Compare target set to reference set. Both either genes or diseases.
  cluster		Compare cluster quality inside clustering. Either genes or diseases.
supported background models
  complete		Random ids will be picked fully randomized.
  term-pres		Random ids will preserve the number of mapped terms for the replaced ids.
 ```
### Result
The validation returns the complete result in a json file
```python
{'status': 'Status text',
 'input_values': {'values': dict(), 'mapped_ids': list()}, 
 'p_values': {'values': dict()}}
```
- **status**: contains either an error message if a mapping failed or "ok" if IDs could be mapped
- **input_values**:
  - **values**: table in dict format with the functional or genetic relevance score(s) determined for solely their input
  - **mapped_ids**: list containing the IDs with non empty annotations per functional or genetic annotation type
- **p_values**: table in dict format with the calculated empirical P-values using the selected background model and other parameters that indicate the significance of the calculated relevance scores derived from the input

As well as separate table files in .csv format for **p_values** and the **relevance score(s)** saved in **input_values**.

If you set the flag `-p` you will also get plots for each type in **p_value** and a 
visualization of the mappability information saved under **mapped_ids**.
### Run with python package
We also offer a [python package](https://pypi.org/project/biodigest).
Check out the [tutorial](https://github.com/bionetslab/digest-tutorial) to see examples of usage in a script.





# Installation

Install conda environment as follows (there also exists a environment.yml but it contains more packages than necessary)
```bash
conda create --name biosteiner python=3.7
conda activate biosteiner
conda install numpy matplotlib pandas networkx pip jupyter
pip install pcst_fast
```

# Running ROBUST

You can simply run robust by calling
```bash
python robust.py data/human_annotated_PPIs_brain.txt data/ms_seeds.txt ms.graphml 0.25 0.9 30 0.1
```
The positional arguments are:
```
[1] file providing the network in the form of an edgelist 
    (tab-separated table, columns 1 & 2 will be used)
[2] file with the seed genes (if table contains more than 
    one column they must be tab-separated; the first column 
    will be used only)
[3] path to output file
[4] initial fraction
[5] reduction factor
[6] number of steiner trees to be computed
[7] threshold
```

The suffix of the path to the output file you specify, determine the format of the output.
You can either choose
- .graphml: A .graphml file is written that contains the following vertex properties: isSeed, significance, nrOfOccurrences, connected_components_id, trees
- .csv: A .csv file which contains a vertex table with #occurrences, %occurrences, terminal (isSeed) 
- everything else: An edge list

# Evaluating ROBUST

For a large-scale empirical evaluation of ROBUST, please follow the instructions given here: https://github.com/bionetslab/robust-eval.

# Citing ROBUST

Please cite ROBUST as follows:
- J. Bernett, D. Krupke, S. Sadegh1, J. Baumbach, S. P. Fekete, T. Kacprowski, M. List1, D. B. Blumenthal: Robust disease module mining via enumeration of diverse prize-collecting Steiner trees, *Bioinformatics* 38(6), pp. 1600-1606, https://doi.org/10.1093/bioinformatics/btab876.
