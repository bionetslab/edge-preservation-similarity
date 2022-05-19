# Edge-Preservation-Similarity
This git provides an exact and an approximated algorithm for computing the edge-preservation similarity between rooted, unordered, node-labeled trees. 

In order to be able to compute the edge-preservation-similarity or run any test first gurobi has to be installed.

## GUROBI

1. register for gurobi
```
https://pages.gurobi.com/registration
```
2. Get free academic licence
```
link to explanation https://www.gurobi.com/academia/academic-program-and-licenses/
```
2.1 Requesting an academic licence
```
https://www.gurobi.com/downloads/end-user-license-agreement-academic/
```
2.2 Accept the conditions

2.3 Follow instuctions at bottom of page and put the grbgetkey with the given numbers in your terminal to install the licence

## Edge-preservation-similarity usage with terminal

For usage of the edge-preservation-similarity with a CLI:

1. Install git
```
pip install git
```
2. Clone this repository
```
Note: This doesn't work with anonymous GitHub. You can either usea workaround(https://github.com/ShoufaChen/clone-anonymous4open) or download the files manually. If you do so please make sure you keep the folder structure as is.
```
3. Install requirements
```
pip install -r requirements.txt
```
4. Go to folder edge_preservation_similarity

5. Run in terminal
```
usage: python CLI_eps.py [required arguments] [optional arguments]

required arguments:
        path for output     string containing the path for the location of the output
        path to trees       string containing one of two options:  
                              - path to .txt file with paths to .gml tree files in each line
                              - paths to all .gml tree files in a row
optional arguments: 
        --algorithm         possibility to choose version of algorithm, choices: {approx,exact}
                                (default: approx)
        --time_limit        possibility to set time limit in seconds for exact algorithm
                                (default: 0 meaning no time limit), data type: int
        --normalize         flag to normalize similarity by dividing by max nr. of edges in tree1 and tree2
                                (default: false, meaning no normalization)
        --both_directions   flag to compute similarity between trees in both directions for more    
                                precise output (default: false, meaning just one direction)
        -h, --help          show this help message and exit
```
### Result
The results are given in two .csv files:
```
similarity_name_of_algorithm.csv    containing the edge-preservation-similarity values
duration_name_of_algorithm.csv      containing the durations for every computation
```



## Compute all tests usage with terminal

The tests are split in 2 parts, scalability tests and validation tests.

### Scalability tests

This test should provide information about the behaviour of the runtime of the edge-preservation-similarity in terms of the tree size. Using the default options of this test the results of the paper will be reproduced.

For usage of the scalability tests of the edge-preservation-similarity with a CLI:

1. Install git
```
pip install git
```
2. Clone this repository
```
Note: This doesn't work with anonymous GitHub. You can either usea workaround(https://github.com/ShoufaChen/clone-anonymous4open) or download the files manually. If you do so please make sure you keep the folder structure as is.
```
3. Install requirements
```
pip install -r requirements.txt
```
4. Go to folder tests

5. Run in terminal
```
usage: scalability_test.py [required arguments] [optional arguments]

CLI for the scalability tests of the edge-preservation-similarity

required arguments:
        output_path         Path to folder where output should be saved

optional arguments:
        --new_trees         flag to compute new scalability trees (default: false, meaning results from paper will be reproduced)
        --time_limit LIMIT  Set time limit in seconds for exact eps algorithm (default: 0 meaning no time limit)
        -h, --help          show this help message and exit

```

### Validation tests

This test compares the edge-preservation-similarity (exact and approximated version) to the tree-edit-distance. As data RNA secondary structure trees are used as described in Wang et al. (2020). First, the similarities between the trees are computed and normalized. Second, the normalized similarities are compared to the functional similarities (Jaccard) of provided Gene Ontology terms.

1. Install git
```
pip install git
```
2. Clone this repository
```
Note: This doesn't work with anonymous GitHub. You can either usea workaround(https://github.com/ShoufaChen/clone-anonymous4open) or download the files manually. If you do so please make sure you keep the folder structure as is.
```
3. Install requirements
```
pip install -r requirements.txt
```
4. Go to folder tests

5. Run in terminal
```
usage: validation_test.py [required arguments] [optional arguments]

CLI for the validation tests of the edge-preservation-similarity

required arguments:
        output_path         Path to folder where output should be saved

optional arguments:
        -h, --help          show this help message and exit
```