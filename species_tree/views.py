from django.shortcuts import render,render_to_response
from .models import Records
from django.http import HttpResponse,HttpResponseRedirect
#from django.template import loader
from datetime import datetime

import os,sys
from Bio.Align.Applications import ClustalwCommandline
#from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.Applications import PhymlCommandline
#from Bio import AlignIO
from Bio import Phylo
import itertools,pydot
#from rpy2.robjects import DataFrame,Matrix
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
from Bio import AlignIO
import copy
ape = importr('ape')
base = importr('base')

def display(request):
    return render(request,"display.html")
def home_page(request):
    return render(request,"home.html")

def save_post(request):

    def juge_os_and_set_PATH():
        if os.name == 'posix' and sys.version_info[0] == 3:
            os.environ['PATH'] = os.environ['PATH'] + ':' + os.getcwd() + '/sowftwares'
        else:
            os.environ['PATH'] = os.environ['PATH'] + ';' + os.getcwd() + '/sowftwares'

    def clustal2phy(file_name_with_path):
        AlignIO.convert(file_name_with_path + '.aln', "clustal", file_name_with_path + ".phy", "phylip")
        return 0

    def construc_tree(file_name_with_path,file_name):
        phyml = PhymlCommandline(input=file_name_with_path + '.phy',outfile=file_name+".nwk")
        phyml()
        return 0

    def parse_alignment(file_name_with_path):
        alignment = AlignIO.read(file_name_with_path + '.aln', 'clustal')
        align_dataframe = pd.DataFrame(data=[list(rec) for rec in alignment], index=[i.id for i in alignment._records],
                                       dtype=np.character)
        return align_dataframe

    def compute_pairwise_distance(file_name_with_path, model='K80'):
        r_alignments_data = ape.read_dna(file_name_with_path + '.aln', format="clustal")
        '''
        model   a character string specifying the evolutionary model to be used;
                must be one of "raw", "N", "TS", "TV", "JC69", "K80" (the default),
                "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock".
        '''
        distance = ape.dist_dna_(r_alignments_data, model=model, as_matrix=True)
        index = list(pandas2ri.ri2py_listvector(distance.rownames))
        distance_numpy = pandas2ri.ri2py_listvector(distance)
        distance_dataframe = pd.DataFrame(distance_numpy, index=index, columns=index)
        return distance_dataframe

    def parse_tree(file_name_with_path,distance_dataframe):
        tree = Phylo.read(file_name_with_path + ".phy_phyml_stats.txt", "newick")
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

    def plot(results):
        # 概率分布直方图
        x = results
        n, bins, patches = plt.hist(x, bins=40, normed=1, histtype='bar', facecolor='green', alpha=0.75)
        plt.title(r'frequency distribution histogram of distance')
        plt.show()
        return 0

    def modify_tree(file_name_with_path,file_name,distance_dataframe, min_number):
        tree = Phylo.read(file_name_with_path + ".phy_phyml_stats.txt", "newick")
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
        Phylo.draw_ascii(newtree)
        Phylo.write(newtree, file_name + '_modify_tree.nwk', 'newick')
        return 0

    def generate_tree(record):

        juge_os_and_set_PATH()
        infile_path = record.inputfile.path
        file_name = record.inputfile.name.split('.')[0]
        file_name_with_path = infile_path.split('.')[0]
        cline = ClustalwCommandline("clustalw2",infile= infile_path,outfile = file_name + ".aln") #Alignment multisequence
        #muscle_cline = MuscleCommandline(input=infile_name + ".fasta", out=infile_name + ".aln",
        #                                 clwstrict=True)  # Alignment multisequence
        cline()

        clustal2phy(file_name_with_path)
        construc_tree(file_name_with_path,file_name)

        distance_dataframe = compute_pairwise_distance(file_name_with_path)
        results = parse_tree(distance_dataframe)
        plot(results)

        modify_tree(file_name_with_path,file_name,distance_dataframe, 0.025)



    if request.method == "POST":
        #try:
        user_name = request.POST.get("user")
        inputfile = request.FILES["inputfile"]
        email = request.POST.get("email")
        record = Records.objects.create(user=user_name, inputfile=inputfile, submit_date=datetime.now(), email=email)
        # print("Submit successfully,waiting for minites to redirect result page.")
        generate_tree(record)
        return HttpResponseRedirect('/species_tree/display/')
        #except:
    # return HttpResponse("Please input completely the information,try again!")
    else:
        return HttpResponse("Input error,please again!")
