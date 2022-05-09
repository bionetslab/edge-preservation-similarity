# Edge-Preservation-Similarity
This git provides an exact and an approximated algorithm for computing the edge-preservation similarity between rooted, unordered, node-labeled trees.

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
git clone https://github.com/bionetslab/edge-preservation-similarity.git
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
        --time_limit=       possibility to set time limit in seconds for exact algorithm (default: 0 meaning no time
                            limit), data type: int
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
TODO