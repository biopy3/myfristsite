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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import AlignIO
import copy,math,pyRserve,time
import zipfile

def every_file_complete_path(dir_path):
    li = []
    for dirpath,dirnames,filenames in os.walk(dir_path):
        for filename in filenames:
            postfix = os.path.splitext(filename)[1]
            if postfix == '.fasta' or postfix == '.fas':
                li.append(os.path.join(dirpath,filename))
    return li

def handle_file(file_name_with_path,infile_path):
    dir_path = file_name_with_path + '_' + time.strftime('%Y%m%d-%H-%M')
    postfix = os.path.splitext(infile_path)[-1]
    file_name = infile_path.split('/')[-1]
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    if postfix == '.zip' :
        zip_file = zipfile.ZipFile(infile_path)
        for names in zip_file.namelist():
            zip_file.extract(names,dir_path)

    if postfix == '.fasta' or postfix == '.fas':
        os.rename(infile_path,dir_path + '/' + file_name)

    return dir_path

def juge_os_and_set_PATH():
    if os.name == 'posix' and sys.version_info[0] == 3:
        os.environ['PATH'] = os.environ['PATH'] + ':' + os.getcwd() + '/species_tree/sowftwares'
    else:
        os.environ['PATH'] = os.environ['PATH'] + ';' + os.getcwd() + '/species_tree/sowftwares'

def clustal2phy(file_name_with_path):
    align = AlignIO.read(file_name_with_path + ".aln","clustal")
    length = 2
    dict = {}
    for record in align:
        if len(record.id) > 10:
            id_list = record.id.split("_")
            id = ""
            for i in range(length):
                id = id + id_list[i][0]
            id = id + id_list[-1][-8:]
            dict[id] = record.id
            record.id = id
    AlignIO.write(align,file_name_with_path + ".phy","phylip")

    return dict

def construc_tree(file_name_with_path, file_name,dict):
    phyml = PhymlCommandline(input=file_name_with_path + '.phy')
    phyml()
    tree = Phylo.read(file_name_with_path + ".phy_phyml_tree.txt", "newick")
    for leaf in tree.get_terminals():
        leaf.name = dict[leaf.name]
    Phylo.write(tree, file_name_with_path + ".phy_phyml_tree.txt", "newick")
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

def plot(results,file_name_with_path):
    # 概率分布直方图
    x = results
    bins = math.ceil(max(results)/0.005)
    n,bins,patches = plt.hist(x, bins=bins, normed=1, histtype='bar', facecolor='green', alpha=0.75)
    plt.title("Frequency distribution of K2P genetic distances \n obtained from successive sister-clade pairwise.")
    plt.savefig(file_name_with_path+'.png',format='png')
    for i in range(len(n)):
        if n[i] == 0:
            break
    return (bins[i]+bins[i+1])/2

def modify_tree(file_name_with_path, file_name, distance_dataframe, min_number):
    tree = Phylo.read(file_name_with_path + ".phy_phyml_tree.txt", "newick")
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

def list_spcies(file_name_with_path):
    f = open(file_name_with_path + "_species_list.csv",'w+',encoding='utf-8')
    tree  = Phylo.read(file_name_with_path + "_modified_tree.nwk",'newick')
    clades = tree.get_nonterminals()
    for clade in clades:
        if clade.is_preterminal:
            li = clade.get_terminals()
            for leaf in li[:-1]:
                f.write(leaf.name + ',')
            f.write(li[-1].name + '\n')
    f.close()
    return 0

@shared_task
def generate_tree(infile_path,send_email,user_name):
    file_name_with_path = os.path.splitext(infile_path)[0]
    dir_path = handle_file(file_name_with_path,infile_path)
    juge_os_and_set_PATH()
    file_path_list = every_file_complete_path(dir_path)
    for i in file_path_list:
        infile_path = i
        file_name_with_path = os.path.splitext(i)[0]
        file_name = '/species_tree/recordsfile/' + file_name_with_path.split('/')[-2] + file_name_with_path.split('/')[-1]
            
        cline = ClustalwCommandline("clustalw2", infile=infile_path,
                                outfile=file_name + ".aln")  # Alignment multisequence
        cline()

        dict = clustal2phy(file_name_with_path)
        construc_tree(file_name_with_path, file_name,dict)


        conn = pyRserve.connect(host='localhost', port=6311)

        distance_dataframe = compute_pairwise_distance(conn,file_name_with_path,'k80')

        matrix_path = file_name_with_path + '_distance_matrix.csv'
        distance_dataframe.to_csv(matrix_path)

        results = parse_tree(file_name_with_path, distance_dataframe)

        divide_line = plot(results, file_name_with_path)

        modify_tree(file_name_with_path, file_name, distance_dataframe, divide_line)

        list_spcies(file_name_with_path)

        # send email
        from_email = settings.DEFAULT_FROM_EMAIL
        email = EmailMessage(
            subject='Hello,' + user_name + ':',
            body='Thank you use the SCPC web service,we send this email with results for you.',
            from_email=from_email,
            to=[send_email]
        )

        email.attach_file('' + file_name + '_modified_tree.nwk')
        email.attach_file('' + file_name + '.png')
        email.attach_file('' + file_name + '.phy_phyml_tree.txt')
        email.attach_file('' + file_name + '_species_list.csv')
        email.attach_file('' + file_name + '_distance_matrix.csv')
        email.send()
        conn.close()

        return 0

