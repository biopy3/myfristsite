# Create your tasks here
from __future__ import absolute_import, unicode_literals
from celery import shared_task
from .models import Records,clustalx_model
from django.core.mail import EmailMessage
from django.conf import settings
import os,sys
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
import zipfile,shutil
import subprocess,os

def every_file_complete_path(dir_path):
    li = []
    try:
        shutil.rmtree(dir_path + '/__MACOSX')
    except:
        pass
    for dirpath,dirnames,filenames in os.walk(dir_path):
        for filename in filenames:
            postfix = os.path.splitext(filename)[1]
            if postfix == '.fasta' or postfix == '.fas' or postfix =='.aln':
                li.append(os.path.join(dirpath,filename))
    return li

def handle_file(file_name_with_path,infile_path):
    dir_path = file_name_with_path + '_' + time.strftime('%Y%m%d-%H-%M')
    postfix = os.path.splitext(infile_path)[-1]
    file_name = infile_path.split('/')[-1]
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    if zipfile.is_zipfile(infile_path):
        zip_file = zipfile.ZipFile(infile_path)
        for names in zip_file.namelist():
            zip_file.extract(names,dir_path)

    if postfix == '.fasta' or postfix == '.fas' or postfix == '.aln':
        os.rename(infile_path,dir_path + '/' + file_name) #move the file

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
            id_ = record.id[-10:]
            dict[id_] = record.id
            record.id = id_
        else:
            dict[record.id] = record.id
    AlignIO.write(align,file_name_with_path + ".phy","phylip")

    return dict

def fasta2phy(file_name_with_path,postfix):
    #try:
    align = AlignIO.read(file_name_with_path + postfix,"fasta")
    #except:
    #  align = AlignIO.read(file_name_with_path + ".fasta","fasta")
    length = 2
    dict = {}
    for record in align:
        if len(record.id) > 10:
            id_ = record.id[-10:]
            dict[id_] = record.id
            record.id = id_
        else:
            dict[record.id] = record.id
    AlignIO.write(align,file_name_with_path + ".phy","phylip")

    return dict

def construc_tree(file_name_with_path,file_name,dict):
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

def compute_pairwise_distance(conn,file_name_with_path,model='K80',postfix='.aln'):
    conn.r.file = file_name_with_path + postfix

#        model   a character string specifying the evolutionary model to be used;
 #               must be one of "raw", "N", "TS", "TV", "JC69", "K80" (the default),
  #              "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock".
    conn.r.model = model
    r_script = '''
        library('ape')
        compute_distance <- function(file,model)
        {
            alignment_data <- read.dna(file=file,format='fasta')
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
    clades = tree.get_nonterminals()
    start_leaives = []
    for clade in clades:
        if clade.is_preterminal():
            start_leaives.append(clade)
    sisters_distance = []
    for clade in start_leaives:
        while tree.root.get_path(clade):
            if len(tree.root.get_path(clade)) >= 2:
                pairs = list(itertools.combinations(tree.root.get_path(clade)[-2].clades, 2))
                for pair in pairs:
                    pair_distance = []
                    for leaf_x in pair[0].get_terminals():
                        leaf_x_name = leaf_x.name
                        for leaf_y in pair[1].get_terminals():
                            leaf_y_name = leaf_y.name
                            pair_distance.append(distance_dataframe[leaf_x_name][leaf_y_name])
                    sisters_distance.append(sum(pair_distance) / len(pair_distance))
                clade = tree.root.get_path(clade)[-2]
            else:  # Only case ::len(newtree.root.get_path(clade)) == 1
                pairs = list(itertools.combinations(tree.root.clades, 2))
                for pair in pairs:
                    pair_distance = []
                    for leaf_x in pair[0].get_terminals():
                        leaf_x_name = leaf_x.name
                        for leaf_y in pair[1].get_terminals():
                            leaf_y_name = leaf_y.name
                            pair_distance.append(distance_dataframe[leaf_x_name][leaf_y_name])
                    sisters_distance.append(sum(pair_distance) / len(pair_distance))
                clade = tree.root
    return sisters_distance

def plot(results,file_name_with_path):
    # 概率分布直方图
    x = results
    bins = math.ceil(max(results)/0.005)
    n,bins,patches = plt.hist(x, bins=bins, normed=1, histtype='bar', facecolor='blue', alpha=1)
    plt.title("Frequency distribution of K2P genetic distances \n obtained from successive sister-clade pairwise.")
    plt.xlabel("genetic distance")
    plt.ylabel("frequency")
    plt.savefig(file_name_with_path+'.png',format='png')
    plt.close()
    
    for i in range(len(n)):
        if n[i] == 0 and n[i+1] != 0:
            break
    return (bins[i+1])

def modify_tree(file_name_with_path, file_name, distance_dataframe, min_number):
    tree = Phylo.read(file_name_with_path + ".phy_phyml_tree.txt", "newick")
    newtree = copy.deepcopy(tree)
    clades = newtree.get_nonterminals()
    start_leaives = []
    for clade in clades:
        if clade.is_preterminal():
            start_leaives.append(clade)
    for clade in start_leaives:
        while newtree.root.get_path(clade):
            if len(newtree.root.get_path(clade)) >= 2:
                pairs = list(itertools.combinations(newtree.root.get_path(clade)[-2].clades, 2))
                distance = []
                for pair in pairs:
                    pair_distance = []
                    for leaf_x in pair[0].get_terminals():
                        leaf_x_name = leaf_x.name
                        for leaf_y in pair[1].get_terminals():
                            leaf_y_name = leaf_y.name
                            pair_distance.append(distance_dataframe[leaf_x_name][leaf_y_name])
                    distance.append(sum(pair_distance) / len(pair_distance))
                if sum(distance) / len(distance) > min_number and len(distance) > 0:
                    for child in newtree.root.get_path(clade)[-2].clades:
                        child.collapse_all()
                    break
                else:
                    clade = newtree.root.get_path(clade)[-2]
            else:  #Only case ::len(newtree.root.get_path(clade)) == 1
                pairs = list(itertools.combinations(newtree.root.clades, 2))
                distance = []
                for pair in pairs:
                    pair_distance = []
                    for leaf_x in pair[0].get_terminals():
                        leaf_x_name = leaf_x.name
                        for leaf_y in pair[1].get_terminals():
                            leaf_y_name = leaf_y.name
                            pair_distance.append(distance_dataframe[leaf_x_name][leaf_y_name])
                    distance.append(sum(pair_distance) / len(pair_distance))
                if sum(distance) / len(distance) > min_number and len(distance) > 0:
                    for child in newtree.root.clades:
                        child.collapse_all()
                    break
                else:
                    clade = newtree.root
    Phylo.write(newtree, file_name_with_path + '_modified_tree.nwk', 'newick')
    return 0

def list_spcies(file_name_with_path):
    f = open(file_name_with_path + "_species_list.csv",'w+',encoding='utf-8')
    tree  = Phylo.read(file_name_with_path + "_modified_tree.nwk",'newick')
    clades = tree.get_nonterminals()
    start_leaives = []
    for clade in clades:
        if clade.is_preterminal():
            start_leaives.append(clade)
    for clade in start_leaives:
        while tree.root.get_path(clade):
            if len(tree.root.get_path(clade)) >= 2:
                for i in tree.root.get_path(clade)[-2].clades.remove(clade):
                    if i.is_terminal():
                        f.write(i.name + '\n')
                clade = tree.root.get_path(clade)[-2]
            else:
                for i in tree.root.clades:
                    if i.is_terminal():
                        f.write(i.name + '\n')
                clade = tree.root
    f.close()


    return 0

def plot_divide_line(list,dir_path):
    plt.scatter(list,list,c='b',marker = 'o')
    plt.savefig(dir_path+'/scatter.png',format='png')
    plt.close()

@shared_task
def generate_tree(infile_path,send_email,user_name,access_code,model):
    file_name_with_path = os.path.splitext(infile_path)[0]
    dir_path = handle_file(file_name_with_path,infile_path)
    juge_os_and_set_PATH()
    file_path_list = every_file_complete_path(dir_path)
    divide_line_list = []
    error_file = []
    successed_file = []
    for i in file_path_list:
        try:
            infile_path = i
            file_name_with_path = os.path.splitext(i)[0]

            file_name = 'species_tree/recordsfile/' + file_name_with_path.split('species_tree/recordsfile/')[-1]
            postfix = os.path.splitext(i)[1]
            if  postfix == '.fasta' or postfix == '.fas':
                dict = fasta2phy(file_name_with_path,postfix)
            else:
                dict = clustal2phy(file_name_with_path)
            construc_tree(file_name_with_path, file_name,dict)

            conn = pyRserve.connect(host='localhost', port=6311)

            distance_dataframe = compute_pairwise_distance(conn,file_name_with_path,model,postfix)

            matrix_path = file_name_with_path + '_distance_matrix.csv'
            distance_dataframe.to_csv(matrix_path)

            results = parse_tree(file_name_with_path, distance_dataframe)
            #save bar_data
            bar_data_path = file_name_with_path+"_bar_data.csv"
            bar_data = np.array(results)
            bar_data.shape = (bar_data.shape[0], 1)
            np.savetxt(bar_data_path, bar_data, delimiter=',')

            divide_line = plot(results, file_name_with_path)
            divide_line_list.append(divide_line)
            modify_tree(file_name_with_path, file_name, distance_dataframe, divide_line)

            list_spcies(file_name_with_path)
            successed_file.append(i)
        except:
            error_file.append(i)
            continue

    plot_divide_line(divide_line_list,dir_path)
    f = open(dir_path+"/error_file.txt",'w')
    for i in error_file:
        f.write(i+'\n')
    f.close()
    f = open(dir_path+"/successed_file.txt",'w')
    for i in successed_file:
        f.write(i+'\n')
    f.close()
    shutil.make_archive(dir_path,'zip',dir_path)
    record = Records.objects.get(access_code=access_code)
    record.resultfile = dir_path + '.zip'
    record.save()
    
    # send email
    from_email = settings.DEFAULT_FROM_EMAIL
    email = EmailMessage(
        subject='Hello,' + user_name + ':',
        body='<p>Thank you use the SCPC web service,we send this email with results for you.Please visit the url:\n</p>\
        </br><a href=http://45.76.122.117:8000/home/result/download/'+access_code+'>\
        http://45.76.122.117:8000/home/result/download/' + access_code + '</a>',
        from_email=from_email,
        to=[send_email]
    )
    email.content_subtype = "html"  # Main content is now text/html
    email.send()
    conn.close()

    return 0

@shared_task
def clustalx_for_align(infile_path,send_email,user_name,access_code):
    file_name_with_path = os.path.splitext(infile_path)[0]
    dir_path = handle_file(file_name_with_path,infile_path)
    file_path_list = every_file_complete_path(dir_path)
    juge_os_and_set_PATH()

    for file_name in file_path_list:
        subprocess.call(["clustalw2","-INFILE="+file_name,"-ALIGN","-QUIET","-OUTPUT=FASTA","-OUTFILE="+ os.path.splitext(file_name)[0]+'_aligned.fasta'])
        os.remove(os.path.splitext(file_name)[0]+'.dnd')
        os.remove(file_name)


    shutil.make_archive(dir_path,'zip',dir_path)
    shutil.rmtree(dir_path)
    record = clustalx_model.objects.get(access_code=access_code)
    record.output_file = dir_path + '.zip'
    record.save()

    #send email
    from_email = settings.DEFAULT_FROM_EMAIL
    email = EmailMessage(
        subject='Hello,' + user_name + ':',
        body='<p>Thank you use the SCPC web service,we send this email with results for you.Please visit the url:\n</p>\
        </br><a href=http://45.76.122.117:8000/clustalx/result_download/'+access_code+'>\
        http://45.76.122.117:8000/clustalx/result_download/' + access_code + '</a>',
        from_email=from_email,
        to=[send_email]
    )
    email.content_subtype = "html"  # Main content is now text/html
    email.send()
