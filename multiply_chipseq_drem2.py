import glob
from pybedtools import BedTool
import unittest
from itertools import chain
import pandas as pd
import numpy as np
import re


class TestStringMethods(unittest.TestCase):

    def test_combine(self):
        list = [("a", [1, 2, 3, 4, 5]),
                ("a", [1, 2, 3, 4, 5, 6]),
                ("a", [1, 2, 3, 4, 5, 6, 7]),
                ("a", [1, 2, 3, 4, 5, 6, 7, 8]),
                ("b", [3, 4, 5, 6, 7, 8, 9]),
                ("b", [1, 2, 3, 4, 5]),
                ("c", [1, 2, 3, 4, 5]),
                ("c", [1, 2, 3, 4, 10]),
                ("c", [1, 2, 3, 4, 11])]
        p_list1={
            "a" : {1, 2, 3, 4, 5, 6, 7, 8},
            "b" : {1, 2, 3, 4, 5, 6, 7, 8, 9},
            "c" : {1, 2, 3, 4, 5, 10, 11}
        }
        self.assertEqual(combine(list),p_list1)


def experi_parser(dir):
    directory=dir + "*.gff3.gz"
    file_list=glob.glob(directory)
    context_list=[]
    tf_list=[]
    i=0
    for file_name in file_list:
        file_name=file_name.split("/")[-1]
        complex=file_name.split(":")
        if len(complex)==8:
            tf, condition, experi, rep, non, ws_version, id, file_type = complex
        else:
            tf, condition, experi, rep, non, ws_version, id= complex
        # covert proper name
        if "-" in tf:      
            match = re.search("\-[0-9]+[A-Za-z]$",tf)  #### eg. LIN-15B --lin-15B
            if match:
                tf=tf[:-1].lower()+tf[-1]
            else:
                tf=tf.lower()
        
        D_Stage, Strain,Temperature=condition.split("#")
        d_Stage = D_Stage.split("=")[1]
        strain = Strain.split("=")[1]
        temperature = Temperature.split("=")[1]
        context=[file_name,rep,d_Stage,strain,temperature]
        context_list.append(context)
        tf_list.append(tf)
    tf_dic = zip(tf_list,context_list)
    return(tf_dic)


def find_target_gene(tf_dic):

    ws220 = '/home/yuxh/zheDA/ChIP-seq/c_elegans_WS220_annotations.bed'
    # ws220 is from c_elegans_WS220_annotations.gff2, exacting gene and pesudogene
    dir="/home/yuxh/zheDA/ChIP-seq/TF/"
    tf_list=[]
    tf_taget_list=[]
    tf_list=[]
    for tf,context in tf_dic:
        tf_chip=dir+context[0]
        a = BedTool(ws220)
        #print(a[1].name)
        result=a.window(tf_chip,u=True,l=1000,r=0)
        #result.head()
        #interval=result[1] # get all bed interval
        gene_list=[iterval.name for iterval in result]
        tf_list.append(tf)
        tf_taget_list.append(gene_list)
    tf_target=list(zip(tf_list,tf_taget_list))
    tf_target.sort()
    return(tf_target)


def combine(tf_target_dic):
    tf_tmp = ""
    tf_target_tmp=[]
    i = 0
    new_list=[]
    for tf,tf_taget_list in tf_target_dic:
        if i == 0:
            i += 1
            tf_tmp = tf
            tf_target_tmp = set(tf_taget_list)

        elif tf_tmp == tf:
            i += 1
            # tf_target_tmp=set(tf_target_tmp).union(set(tf_taget_list))
            # this can be union or intersection
            tf_target_tmp=set(tf_target_tmp).intersection(set(tf_taget_list))
            if i==len(tf_target_dic):
                new_list.append((tf_tmp, tf_target_tmp))
        else:
            i += 1
            new_list.append((tf_tmp,tf_target_tmp))
            tf_tmp = tf
            tf_target_tmp = set(tf_taget_list)

            if i==len(tf_target_dic):
                new_list.append((tf_tmp, tf_target_tmp))
    return(dict(new_list))

def generate_matrix(dic):

    # get all target gene
    # all_target=set(list(chain(*([target for target in dic.values()]))))
    # replace all_target_gene=dic.values() with all possible gene_id
    gene_id=open("/home/yuxh/zheDA/ChIP-seq/c_elegans.WS220.geneIDs.txt","r")
    all_target=[]
    for lines in gene_id:
        line=lines.strip("\n").split(",")
        if "" in line:
            line.remove("")
        gene=";".join(line)
        all_target.append(gene)

    bool_matrix=[]
    tf_matrix=[]
    for tf, target in dic.items():
        for gene_lt in all_target:
            list_tmp=gene_lt.split(";")
            # any gene in target
            if any([(x in target) for x in list_tmp]):
                bool_matrix.append(str(1))
            else:
                bool_matrix.append(str(0))
        tf_matrix.append((tf,bool_matrix))
        bool_matrix=[]

    new_matrix=pd.DataFrame(data=dict(tf_matrix),index=all_target)
    new_matrix.to_csv("c.elegans_tf_union_drem2.csv",sep="\t")



if __name__ == '__main__':
    #unittest.main()

    tf_dic=experi_parser("/home/yuxh/zheDA/ChIP-seq/TF/")
    tf_target=find_target_gene(tf_dic)
    new_list=combine(tf_target)
    generate_matrix(new_list)





