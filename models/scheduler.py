# -*- coding: utf-8 -*-
#sched
import networkx as nx
from networkx.readwrite import json_graph
import json
import fishertest,fisher_standalone
import itertools,sys,os,math,scipy
import numpy as np
import operator,math
from operator import itemgetter

def load_input_list(file):
    f1=open(file,"rU")
    seq=f1.readline()
    temp=seq.strip().split("\t")
    column=len(temp)
    uniprot_list={}
    mapping={}
    for i in range(column):
        uniprot_list[i]=[]
        mapping[i]=temp[i]
    seq=f1.readline()
    while(seq!=""):
        seq=seq.split("\t")
        for i in range(column):
            if len(seq[i])>1:
                uniprot_list[i].append(seq[i].strip())
        seq=f1.readline()
    for i in uniprot_list:
		uniprot_list[i]=list(set(uniprot_list[i]))
    return uniprot_list,mapping

def expression(folder,tissue):
    
    tissue_expr={}
    tissue_number=float(len(tissue))
    if tissue_number>1:
        for j in tissue:
            f1=open(folder+"/"+j+".txt","r")
            seq=f1.readline()
            
            while(seq!=""):
                seq=seq.strip().split("\t")
                if tissue_expr.has_key(seq[0]):
                    tissue_expr[seq[0]]=tissue_expr[seq[0]]+float(seq[1])
                else:
                    tissue_expr[seq[0]]=[]
                    tissue_expr[seq[0]]=float(seq[1])
                seq=f1.readline()
        
        for i in tissue_expr:
            tissue_expr[i]=tissue_expr[i]/tissue_number
    else:
        f1=open(folder+"/"+tissue[0]+".txt","r")
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            tissue_expr[seq[0]]=float(seq[1])
            seq=f1.readline()
    return tissue_expr

def mcn(nodes,graph_nodes,graph,tissue_expr,graph_choice,coexpression,folder,db_values,annotation):
    path_mean={}
    path_count={}
    seek={}
    path_score={}
    path_selection={}
    nodes_enriched_value={}
    nodes_enriched_path={}
    combination=list(itertools.combinations((list(set(graph_nodes).intersection(set(nodes)))),2))
    for i in nodes:
        seek[i]={}
    for i in combination:
        path_count[i]={}
        path_score[i]=10.0
        path_mean[i]=0.0
        path_selection[i]=[]
        nodes_enriched_value[i]=100
        nodes_enriched_path[i]=[]
        f1=open(folder+"data/path/"+graph_choice+"/index/"+i[0]+".txt","r")
        seq=f1.readline()
        while(seq!=""):
            seq= seq.strip().split("\t")
            if seq[0]==i[1]:
				seek[i[0]][seq[0]]=int(seq[1])
            seq=f1.readline()

    for i in seek:
        if len(seek[i])>0:
            for j in seek[i]:
                f1=open(folder+"data/path/"+graph_choice+"/"+i+".txt")
                f1.seek(seek[i][j])
                seq=f1.readline()
                key_start=seq.split("|")[0][1:].strip()
                key_end=seq.split("|")[1].strip()
                key=(key_start,key_end)
                seq=f1.readline()
                count=0
                while(seq[0]!=">"):
                    path_seq=seq.strip()
                    seq=seq.strip().split()
                    node_list=seq
                    path_exp=[]
                    path_coex=[]
                    path_coex.append(graph[node_list[0]][node_list[1]]["coex"][2])
                    for k in range(len(node_list[1:-1])):
                        path_coex.append(graph[node_list[k+1]][node_list[k+2]]["coex"][2])
                        if tissue_expr.has_key(node_list[k+1]):
                            path_exp.append(tissue_expr[node_list[k+1]])
                        else:
                            path_exp.append(0.0001)
                    seq=f1.readline()
                    if len(path_exp)==0:
						val1=1.00000
                    else:
						val1=scipy.stats.mstats.gmean(path_exp)
                    if coexpression=="yes":
                        val2=scipy.stats.mstats.gmean(path_coex)
                        total_prob=math.sqrt(val1*val2)
                    else:
						total_prob=val1
                    if path_count[key].has_key(total_prob):
						path_count[key][tuple(node_list)]=total_prob
                    else:
						#path_count[key][(node_list)]
						path_count[key][tuple(node_list)]=total_prob
                    
                    path_mean[key]=path_mean[key]+total_prob
                    count=count+1
                path_mean[key]=path_mean[key]/float(count)

    for i in path_count:
        if len(path_count[i])==1:
            nodes_enriched_path[i]=list(path_count[i].keys()[0])
        else:
            
            for k in path_count[i]:
                if path_count[i][k]>= path_mean[i]:

                    temp_C=[]
                    temp_P=[]
                    temp_F=[]
                    temp_R=[]
                    temp_K=[]
                    for j in k:
                        if annotation["P"].has_key(j):
                            temp_P.append(annotation["P"][j])
                        if annotation["K"].has_key(j):
                            temp_K.append(annotation["K"][j])
                        if annotation["R"].has_key(j):
                            temp_R.append(annotation["R"][j])
                        if annotation["C"].has_key(j):
                            temp_C.append(annotation["C"][j])
                        if annotation["F"].has_key(j):
                            temp_F.append(annotation["F"][j])
                    value_p=0.0
                    value_c=0.0
                    value_f=0.0
                    value_r=0.0
                    value_k=0.0
                    temp_P=sum(temp_P,[])
                    temp_C=sum(temp_C,[])
                    temp_F=sum(temp_F,[])
                    temp_K=sum(temp_K,[])
                    temp_R=sum(temp_R,[])
                    total_db_terms=float(len(temp_C)+len(temp_P)+len(temp_R)+len(temp_F)+len(temp_K))
                    temp_P=list(set(temp_P))
                    temp_C=list(set(temp_C))
                    temp_F=list(set(temp_F))
                    temp_K=list(set(temp_K))
                    temp_R=list(set(temp_R))
                    for j in temp_P:
                        if db_values["P"].has_key(j):
                            value_p=value_p+ db_values["P"][j]
                    for j in temp_C:
                        if db_values["C"].has_key(j):
                            value_c=value_c+ db_values["C"][j]
                    for j in temp_R:
                        if db_values["R"].has_key(j):
                            value_r=value_r+ db_values["R"][j]
                    for j in temp_F:
                        if db_values["F"].has_key(j):
                            value_f=value_f+ db_values["F"][j]
                    for j in temp_K:
                        if db_values["K"].has_key(j):
                            value_k=value_k+ db_values["K"][j]

                    annotation_value=(value_p+value_c+value_r+value_f+value_k)/total_db_terms
                    if annotation_value<nodes_enriched_value[i]:
                        nodes_enriched_value[i]=annotation_value
                        nodes_enriched_path[i]=list(k)
   
    return list(set(sum(nodes_enriched_path.values(),[])))

def load_db_values(folder,tissue,folder_annotation):
    db_values={}
    annotation={}
    tissue_number=float(len(tissue))
    for i in ["C","P","F","R","K"]:
        db_values[i]={}
        annotation[i]={}
    for i in tissue:
        for j in ["C","P","F","R","K"]:
            f1=open(folder+"/"+i+"/"+j+"_prob.txt","r")
            seq=f1.readline()
            seq=f1.readline()
            while(seq!=""):
                seq=seq.strip().split("\t")
                if db_values[j].has_key(seq[0]):
                    db_values[j][seq[0]]= db_values[j][seq[0]]+float(seq[1])
                else:
                    db_values[j][seq[0]]=float(seq[1])
                seq=f1.readline()
            f1.close()
    for j in ["C","P","F","R","K"]:
        f1=open(folder_annotation+"/"+j+".txt","r")
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            annotation[j][seq[0]]=seq[1:]
            seq=f1.readline()
        f1.close()
    if len(tissue)>1:
        for i in ["C","P","F","R","K"]:
            for j in db_values[i]:
                db_values[i][j]=db_values[i][j]/tissue_number
    return db_values,annotation


human_proteome_map_adult={0: 'adrenal_gland', 1: 'colon', 2: 'esophagus', 3: 'female_gonad', 4: 'frontal_cortex', 5: 'gall_bladder', 6: 'heart', 7: 'kidney', 8: 'liver', 9: 'lung', 10: 'monocyte',
                          11: 'pancreas', 12: 'prostate_gland', 13: 'rectum', 14: 'retina', 15: 'spinal_cord', 16: 'testis', 17: 'urinary_bladder'}
protein_atlas_cancer={0: 'breast_cancer', 1: 'carcinoid', 2: 'cervical_cancer', 3: 'colorectal_cancer', 4: 'endometrial_cancer', 5: 'glioma', 6: 'head_and_neck_cancer', 7: 'liver_cancer', 8: 'lung_cancer',
                      9: 'lymphoma', 10: 'melanoma', 11: 'ovarian_cancer', 12: 'pancreatic_cancer', 13: 'prostate_cancer', 14: 'renal_cancer', 15: 'skin_cancer', 16: 'stomach_cancer', 17: 'testis_cancer',
                      18: 'thyroid_cancer', 19: 'urothelial_cancer'}
human_proteome_map_fetal={0: 'brain', 1: 'female_gonad', 2: 'gut', 3: 'heart', 4: 'liver', 5: 'placenta'}
peptideatlas={0: 'brain', 1: 'heart', 2: 'kidney', 3: 'liver', 4: 'plasma', 5: 'urine', 6: 'whole_organism'}
protein_atlas_basal={0: 'adrenal_gland', 1: 'appendix', 2: 'bone_marrow', 3: 'breast', 4: 'bronchus', 5: 'cerebellum', 6: 'cerebral_cortex', 7: 'cervix_uterine', 8: 'colon', 9: 'duodenum', 10: 'endometrium_1',
                     11: 'endometrium_2', 12: 'epididymis', 13: 'esophagus', 14: 'fallopian_tube', 15: 'gallbladder', 16: 'heart_muscle', 17: 'hippocampus', 18: 'kidney', 19: 'lateral_ventricle', 20: 'liver',
                     21: 'lung', 22: 'lymph_node', 23: 'nasopharynx', 24: 'oral_mucosa', 25: 'ovary', 26: 'pancreas', 27: 'parathyroid_gland', 28: 'placenta', 29: 'prostate', 30: 'rectum', 31: 'salivary_gland',
                     32: 'seminal_vesicle', 33: 'skeletal_muscle', 34: 'skin_1', 35: 'skin_2', 36: 'small_intestine', 37: 'smooth_muscle', 38: 'soft_tissue_1', 39: 'soft_tissue_2', 40: 'spleen', 41: 'stomach_1',
                     42: 'stomach_2', 43: 'testis', 44: 'thyroid_gland', 45: 'tonsil', 46: 'urinary_bladder', 47: 'vagina'}
proteomicsdb={0: 'ascitic_fluid', 1: 'cardia_of_stomach', 2: 'cerebral_cortex', 3: 'colon', 4: 'earwax', 5: 'esophagus', 6: 'fallopian_tube', 7: 'female_gonad', 8: 'gall_bladder', 9: 'hair_follicle',
              10: 'kidney', 11: 'liver', 12: 'lung', 13: 'lymph_node', 14: 'milk', 15: 'nasopharynx', 16: 'oral_cavity', 17: 'pancreas', 18: 'placenta', 19: 'prostate_gland', 20: 'rectum', 21: 'saliva',
              22: 'saliva_secreting_gland', 23: 'seminal_vesicle', 24: 'skin', 25: 'spleen', 26: 'stomach', 27: 'testis', 28: 'thyroid_gland', 29: 'tonsil', 30: 'uterine_cervix', 31: 'uterus', 32: 'vulva'}
PaxDb_integrated={0: 'brain', 2: 'COLON', 3: 'esophagus', 4: 'female_gonad', 5: 'gall_bladder', 6: 'heart', 7: 'kidney', 8: 'liver', 9: 'lung', 10: 'pancreas', 11: 'placenta',
            12: 'plasma', 13: 'platelet', 14: 'prostate_gland', 15: 'rectum', 16: 'saliva', 17: 'skin', 18: 'testis', 19: 'urine', 20: 'uterus', 21: 'whole_organism'}

def tsfa(uid,file,graph,threshold,db,tissue,coexpression,root_folder,c_w,f_w,p_w,r_w,k_w,background):
    
    expr={"0":"protein_atlas_basal","1":"protein_atlas_cancer","2":"humanproteomemap_adult","3":"humanproteomemap_fetal","4":"PaxDb_integrated","5":"peptideatlas","6":"proteomicsdb"}
    net_db={"0":"intact","1":"intact_high_confidence","2":"biogrid"}
    start_nodes,mapping=load_input_list(file)
    f1=open(root_folder+"static/results/"+uid+"/mapping_column.txt","w")
    common_proteins={}
    for i in mapping:
        f1.write(str(i)+"\t"+mapping[i]+"\n")
    graph_net=net_db[graph]
    nodes_graph={}
    graph=nx.read_gpickle(root_folder+"/data/graph/"+graph_net)
    graph.graph['name']=graph_net
    graph_nodes=graph.nodes()
    expr_db=expr[db]
    db_values,annotation=load_db_values(root_folder+"data/score/"+background+"/"+graph_net+"_db_prob/"+expr_db,tissue,root_folder+"data/annotation/"+background)
    tissue_expr=expression(root_folder+"data/tissue_expression/"+expr_db+"/"+graph_net,tissue)

    f1=open(root_folder+"static/results/"+uid+"/nodes.txt","w")
    f2=open(root_folder+"static/results/"+uid+"/nodes_present_in_graph.txt","w")
    for i in start_nodes:
        
        if not os.path.exists(root_folder+"/static/results/"+uid+"/"+str(i)+"_fisher/"):
            os.makedirs(root_folder+"/static/results/"+uid+"/"+str(i)+"_fisher")
        if not os.path.exists(root_folder+"/static/results/"+uid+"/"+str(i)+"_graph/"):
            os.makedirs(root_folder+"/static/results/"+uid+"/"+str(i)+"_graph")
        f1.write(str(i)+"\t"+"\t".join(start_nodes[i])+"\n")
        
        nodes=list(set(start_nodes[i]).intersection(set(graph_nodes)))
        f2.write(str(i)+"\t"+"\t".join(nodes)+"\n")
		# se il numero di nodi e uguale a uno non procedere
    
        if len(nodes)>1:
           
            #aggiungi i nodi iniziali.
            nodes_graph[i]=mcn(nodes,graph.nodes(),graph,tissue_expr,graph_net,coexpression,root_folder,db_values,annotation)
            graph_mcn=graph.subgraph(nodes_graph[i])
    
    f1=open(root_folder+"static/results/"+uid+"/nodes_graph.txt","w")
    for i in nodes_graph:
        stats=""
        graph.graph['name']=graph_net
        stats=stats+nx.info(graph)+"\n"
        graph_mcn=graph.subgraph(nodes_graph[i])
        graph_mcn.graph['name']="Final Network"
        stats=stats+nx.info(graph_mcn)
        nx.write_gpickle(graph_mcn,root_folder+"static/results/"+uid+"/"+str(i)+"_graph"+"/graph_mcn.gpickle")
        f1.write(str(i)+"\t"+"\t".join(list(set(nodes_graph[i])))+"\n")
        f2=open(root_folder+"/static/results/"+uid+"/"+str(i)+"_graph/not_present.txt",'w')
        f2.write("\t".join(list(set(start_nodes[i]).difference(set(graph_mcn.nodes())))))
        f2.close()
        f2=open(root_folder+"/static/results/"+uid+"/"+str(i)+"_graph/stats.txt",'w')
        f2.write(stats)
        f2.close()
    
def tfa(uid,file,graph_db,threshold,topology_ranking):
    test={}

def order_res(order,bw):
    r={}
    removable=[]
    for i in order:
        if r.has_key(order[i]):
            r[order[i]].append(i)
        else:
            r[order[i]]=[]
            r[order[i]].append(i)
    for i in r:
        removable=removable+sorted(r[i], key=bw.get)
    return removable

def check(path,removable,nodes_list):
	temp_rem={}
	rem_path={}
	path_lenght={}
	rem=[]
	ess=[]
	for j in path:
		temp_rem[j]=[]
		rem_path[j]=[]
		path_lenght[j]=len(path[j])		
	for i in removable:
		for j in path:
			for k in range(path_lenght[j]):
				if i in path[j][k] and k not in rem_path[j]:
				   	temp_rem[j].append(k)
		#check se qualche path e' afflitto
		count=0
		for j in path:
			union=temp_rem[j]+rem_path[j]
			union=list(set(union))
			if path_lenght[j]-len(union)>0:
				count=count+1
		if count==len(path):
			rem.append(i)
			for j in path:
				rem_path[j].extend(temp_rem[j])
				temp_rem[j]=[]
		if count!=len(path):
			ess.append(i)
			for j in path:
				temp_rem[j]=[]
	protein_list=ess+nodes_list
	protein_list=list(set(protein_list))
	return protein_list
