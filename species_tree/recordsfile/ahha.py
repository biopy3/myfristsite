# Create your tasks here
from Bio import Phylo
import itertools
import pandas as pd
from matplotlib import pyplot as plt
import pyRserve

conn = pyRserve.connect(host='localhost', port=6311)
file_name_with_path = '/Users/mac/Documents/myfristsite/species_tree/recordsfile/Callosciurinae'

def compute_pairwise_distance(conn,file_name_with_path, model='K80'):
    conn.r.file = file_name_with_path + '.aln'
    conn.r.model = model
    r_script='''
    library('ape')
    compute_distance <- function(file,model) {
        alignment_data <- read.dna(file=file,format='clustal')
        index <- attr(alignment_data,"dimnames")[[1]]
        distance_data <- dist.dna(alignment_data,model=model,as.matrix=TRUE)
        result <- list(a = distance_data,b = index)
        return(result)
    }
    result <- compute_distance(file,model)
    distance_data <- result$a
    index <- result$b
    '''
    conn.eval(r_script)
    return  conn.r.distance_data,list(conn.r.index)


def plot(results):
    x = results
    plt.hist(x, bins=40, normed=1, histtype='bar', facecolor='green', alpha=0.75)
    plt.title(r'frequency distribution histogram of distance')
    plt.savefig(file_name_with_path + '.png', format='png')
    return 0

def parse_tree(file_name_with_path, distance_dataframe):
    tree = Phylo.read(file_name_with_path + ".phy_phyml_tree.txt", "newick")
    li_0 = []
    for clade in tree.get_nonterminals():
        li = list(itertools.combinations(clade, 2))
        li_1 = []
        for clade_pair in li:
            li_2 = []
            for leaf in clade_pair[0].get_terminals():
                leaf_1 = leaf.name
                for leaf in clade_pair[1].get_terminals():
                    leaf2 = leaf.name
                    li_2.append(distance_dataframe[leaf_1][leaf2])
            li_1.append(sum(li_2) / len(li_2))
        li_0.append(max(li_1))
    return li_0

distance_data , index = compute_pairwise_distance(conn,file_name_with_path)
distance_dataframe = pd.DataFrame(distance_data,
                                      index=list(index), columns=list(index))
results = parse_tree(file_name_with_path, distance_dataframe)
plot(results)

