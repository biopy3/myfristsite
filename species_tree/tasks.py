# Create your tasks here
from __future__ import absolute_import, unicode_literals
from celery import shared_task
from django.core.mail import EmailMessage
from django.conf import settings
import os,sys
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.Applications import PhymlCommandline
from Bio import Phylo
import itertools
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from Bio import AlignIO
import copy
import pyRserve



def juge_os_and_set_PATH():
    if os.name == 'posix' and sys.version_info[0] == 3:
        os.environ['PATH'] = os.environ['PATH'] + ':' + os.getcwd() + '/species_tree/sowftwares'
    else:
        os.environ['PATH'] = os.environ['PATH'] + ';' + os.getcwd() + '/species_tree/sowftwares'

def clustal2phy(file_name_with_path):
    AlignIO.convert(file_name_with_path + '.aln', "clustal", file_name_with_path + ".phy", "phylip")
    return 0

def construc_tree(file_name_with_path, file_name):
    phyml = PhymlCommandline(input=file_name_with_path + '.phy')
    phyml()
    return 0

def parse_alignment(file_name_with_path):
    alignment = AlignIO.read(file_name_with_path + '.aln', 'clustal')
    align_dataframe = pd.DataFrame(data=[list(rec) for rec in alignment], index=[i.id for i in alignment._records],
                                   dtype=np.character)
    return align_dataframe

def compute_pairwise_distance(conn,file_name_with_path, model='K80'):
    conn.r.file = file_name_with_path + '.aln'

#        model   a character string specifying the evolutionary model to be used;
 #               must be one of "raw", "N", "TS", "TV", "JC69", "K80" (the default),
  #              "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock".
    conn.r.model = model
    r_script = '''
        library('ape')
        compute_distance <- function(file,model)
        {
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
    distance_dataframe = pd.DataFrame(conn.r.distance_data,
                                      index=list(conn.r.index), columns=list(conn.r.index))

    return distance_dataframe

def parse_tree(file_name_with_path, distance_dataframe):
    tree = Phylo.read(file_name_with_path + "_phy_phyml_tree.txt", "newick")
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

def plot(results,file_name_with_path):
    # 概率分布直方图
    x = results

    plt.hist(x, bins=40, normed=1, histtype='bar', facecolor='green', alpha=0.75)
    plt.title(r'frequency distribution histogram of distance')
    plt.savefig(file_name_with_path+'.png',format='png')
    return 0

def modify_tree(file_name_with_path, file_name, distance_dataframe, min_number):
    tree = Phylo.read(file_name_with_path + "_phy_phyml_tree.txt", "newick")
    newtree = copy.deepcopy(tree)
    clades = newtree.get_nonterminals()
    for clade in clades[1:]:
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
        if max(li_1) < min_number:
            newtree.collapse(clade)
    Phylo.write(newtree, file_name_with_path + '_modified_tree.nwk', 'newick')
    return 0

@shared_task
def generate_tree(file_name,infile_path):
    juge_os_and_set_PATH()
    file_name_with_path = infile_path.split('.')[0]
    cline = ClustalwCommandline("clustalw2", infile=infile_path,
                                outfile=file_name + ".aln")  # Alignment multisequence
    # muscle_cline = MuscleCommandline(input=infile_name + ".fasta", out=infile_name + ".aln",
    #                                 clwstrict=True)  # Alignment multisequence
    cline()

    clustal2phy(file_name_with_path)
    construc_tree(file_name_with_path, file_name)

    return 0

@shared_task
def modifytree(file_name,infile_path,send_email,user_name):

    conn = pyRserve.connect(host='localhost', port=6311)
    file_name_with_path = infile_path.split('.')[0]

    distance_dataframe = compute_pairwise_distance(conn,file_name_with_path,'k80')

    results = parse_tree(file_name_with_path, distance_dataframe)

    plot(results, file_name_with_path)

    
    modify_tree(file_name_with_path, file_name, distance_dataframe, 0.025)

    # send email
    from_email = settings.DEFAULT_FROM_EMAIL
    email = EmailMessage(
        subject='Hello,' + user_name + ':',
        body='Thank you use the ISDL web service,we send this email with results for you',
        from_email=from_email,
        to=[send_email]
    )

    email.attach_file('' + file_name + '_modified_tree.nwk')
    #email.attach_file('' + file_name + '.png')
    email.send()

    conn.close()

    return 0

