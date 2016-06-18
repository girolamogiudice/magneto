#!/usr/bin/env python
# -*- coding: utf-8 -*-
from gluon import *
import json
from scipy.stats import fisher_exact
import networkx as nx
import sys,os,math
def load(uuid,ids,protein,threshold,comp_list,choice,folder,start_nodes,domain_db,protein_descr,count_annotation,root_second_level,fisher_annotation,sample,sample_number,control,controls):
    bg={"0":"proteome","1":"uniprot"}
    colors={"C":"#d5f0f4","P":"#2171b5","F":"#6baed6","R":"#abe16c","K":"#5f8726","O":"#7d6396",
            "KDr":"#b23131","KDi":"#ffcdff","DB":"#ff0000","Or":"#c19bce","HPi":"#A9A9A9","T":"#b56b19"}
    #drugs http://www.color-hex.com/color-palette/17462
    #disease http://www.color-hex.com/color-palette/17377
    #pathways http://www.color-hex.com/color-palette/16868
    #go terms 
    col={}
    mat={}
    neigh={}
    fisher_json_table=[]
    databases={"C":"Gene Ontology","F":"Gene Ontology","P":"Gene Ontology","R":"Pathway","K":"Pathway","O":"Disease","Or":"Disease","KDi":"Disease","DB":"Drugs","KDr":"Drugs","HPi":"Virus","T":"Toxins"}
    alias_dict={"C":"Component","F":"Function","P":"Process","K":"Kegg Pathway","R":"Reactome","O":"Omim","KDi":"Kegg Disease","KDr":"Kegg Drug","DB":"Drugbank","Or":"Orphanet","HPi":"Virus","T":"Toxins"}
    graph_mcn = nx.read_gpickle(folder+"static/results/"+uuid+"/"+ids+"_graph/graph_mcn.gpickle")
    f1=open(folder+"/static/results/"+uuid+"/"+ids+"_graph/neigh.txt","w")
    for i in graph_mcn.nodes():
        f1.write(i+"\t"+"\t".join(graph_mcn.neighbors(i))+"\n")
        neigh[i]=graph_mcn.neighbors(i)
    f1.close()
    path=folder+"data/annotation/"+choice+"/"
    children={}
    f7=open(folder+"/static/results/"+uuid+"/"+ids+"_graph/cluster_id.txt","w")
    f8=open(folder+"/static/results/"+uuid+"/"+ids+"_bar.txt","w")
    f8.write("db\tStandard\tCommon\tMAGNETO\tTotal\n")
    f9=open(folder+"/static/results/"+uuid+"/"+ids+"_chord_plot.txt","w")
    temp_db={}
    for ii in comp_list:
        f3=open(folder+"/static/results/"+uuid+"/"+ids+"_graph/"+ii+"_word_fisher.txt","w")
        f4=open(folder+"/static/results/"+uuid+"/"+ids+"_graph/"+ii+"_cluster_id.txt","w")
        f5=open(folder+"/static/results/"+uuid+"/"+ids+"_graph/"+ii+"_cluster_db.txt","w")
        f6=open(folder+"/static/results/"+uuid+"/"+ids+"_graph/"+ii+"_descr.txt","w")
        prot_fisher=[]
        temp_db[ii]=[]
        clusterdb={}
        clusterid={}
        
        words_count={}
        f1=open(folder+"data/annotation/"+choice+"_words/"+ii+"_count.txt","r")
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            words_count[seq[0]]=int(seq[1])
            seq=f1.readline()
        f1.close()
        
        words={}
        f1=open(folder+"data/annotation/"+choice+"_words/"+ii+".txt","r")
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            words[seq[0]]=seq[1].split()
            seq=f1.readline()
        descr={}
        f1=open(path+ii+"_descr.txt","r")
        seq=f1.readline()
        while (seq!=""):
            seq=seq.strip().split("\t")
            descr[seq[0]]=seq[1]
            seq=f1.readline()
        f1.close()
        f1=open(path+ii+".txt","r")
        seq=f1.readline()
        fisher={}
        fisher_count={}
        id_uniprot={}
        while(seq!=""):
            seq=seq.strip().split("\t")
            fisher[seq[0]]=seq[1:]
            seq=f1.readline()
        f1.close()
        f1=open(path+ii+"_count.txt","r")
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            fisher_count[seq[0]]=int(seq[1])
            seq=f1.readline()
        f1.close()

        fisherset={}
        fisherset_count=[]
        annotation={}
        annotation_gene={}
        for i in protein:
            if i in fisher:
                for j in fisher[i]:
                    if j in fisherset:
                        annotation_gene[j].append(i+"("+protein_descr[i]+")")
                        annotation[j].append(i)
                        fisherset[j]=fisherset[j]+1
                    else:
                        annotation_gene[j]=[]
                        annotation_gene[j].append(i+"("+protein_descr[i]+")")
                        annotation[j]=[]
                        annotation[j].append(i)
                        fisherset[j]=1
            else:
                fisherset_count.append(i)

        totalfisher=len(fisher)
        numberofproteins=len(protein)-len(list(set(fisherset_count)))
        fisher={}
        fisher_value={}
        nmap={}
        words_count_enriched={}
        count=0
        lenfisherset=len(list(set(fisherset)))
        fisher_json=[]
        
        for i in fisherset:
            a=fisherset[i]
            b=numberofproteins-a
            c=fisher_count[i]-a
            d=totalfisher-a-b-c
            table=[[a,b],[c,d]]
            fisher[i]=fisher_exact(table,alternative ="greater")[1]
            if fisher[i]<(threshold/lenfisherset):
                root_second_level[ii]={"name":databases[ii],"children":[{"name":alias_dict[ii],"children":[]}]}
                if control:
                    controls[ii][i]=""
                if sample_number>1:
                    if count_annotation[ii].has_key(i):
                        count_annotation[ii][i]["size"]=count_annotation[ii][i]["size"]+1
                        count_annotation[ii][i]["sample"]=count_annotation[ii][i]["sample"]+" "+sample
                    else:
                        count_annotation[ii][i]={"name":i,"size":1,"description":descr[i],"proteins":" ".join(annotation_gene[i]),"sample":sample}
                domain_db[ii].append(i)
                nmap[i]={"id":count,"name":i,"description":descr[i]}
                count=count+1
                fisher_value[i]=fisher[i]

        f2=open(folder+"/static/results/"+uuid+"/"+"/"+ids+"_graph/"+ii+"fisher.txt","w")
        f2.write("id\tp_value\tproteins_involved\tdescription\tproteins\n")
        common=0
        magneto=0
        standard=0
        for i in fisher_value:
            if i in fisher_annotation[ii]:
                fisher_json.append({"id":i,"p value":str(fisher_value[i]),"proteins involved":str(len(annotation[i])),"description":descr[i],"proteins":annotation_gene[i],"common":1})
                common=common+1
            else:
                fisher_json.append({"id":i,"p value":str(fisher_value[i]),"proteins involved":str(len(annotation[i])),"description":descr[i],"proteins":annotation_gene[i],"common":0})
                magneto=magneto+1
            f2.write(i+"\t"+str(fisher_value[i])+"\t"+str(len(annotation[i]))+"\t"+descr[i]+"\t"+" ".join(annotation_gene[i])+"\n")
            f6.write(i+"\t"+descr[i]+"\n")
                
                
            for k in words[i]:
                if words_count_enriched.has_key(k):
                    words_count_enriched[k]=words_count_enriched[k]+1
                else:
                    words_count_enriched[k]=1

            for k in annotation[i]:
                prot_fisher.append(k)
                if i in clusterid:
                    clusterid[i].append(k)
                else:
                    clusterid[i]=[]
                    clusterid[i].append(k)
                if k in clusterdb:
                    clusterdb[k].append(i)
                else:
                    nmap[k]={"id":count,"name":k+"("+protein_descr[k]+")","description":protein_descr[k],"descr":protein_descr[k][0]}
                    count=count+1
                    clusterdb[k]=[]
                    clusterdb[k].append(i)


                if k in col:
                    col[k].append(colors[ii])
                else:
                    col[k]=[]
                    col[k].append(colors[ii])
        f2.close()
        for i in fisher_annotation[ii]:
            if not i in fisher_value:
                fisher_json.append({"id":i,"p value":"-1","proteins involved":str(len(annotation[i])),"description":descr[i],"proteins":annotation_gene[i],"common":2})
                standard=standard+1
        if not (standard==0 and common==0 and magneto==0):
            fisher_json_table.append({"database":alias_dict[ii],"Standard":str(standard),"Shared":str(common),"MAGNETO":str(magneto),"Total":str(common+magneto)})
            f8.write(alias_dict[ii]+"\t"+str(standard)+"\t"+str(common)+"\t"+str(magneto)+"\t"+str(common+magneto)+"\n")
        json.dump(fisher_json,open(folder+"/static/results/"+uuid+"/"+"/"+ids+"_graph/"+ii+"fisher.json","w"))
        json.dump(fisher_json_table,open(folder+"/static/results/"+uuid+"/"+ids+"_bar.json","w"))
        number_of_word=len(words_count_enriched)
        total_word=len(words_count)
        word_fisher=[]
        for i in words_count_enriched:
            a=words_count_enriched[i]
            #############################################
            ##DEvo fare un recheck di questa funzione####
            #############################################
            #print i,words_count_enriched[i]
            b=words_count_enriched[i]-a
            c=words_count[i]-a
            d=sum(words_count.values())-a-b-c
            fisher[i]=fisher_exact(table,alternative ="greater")[1]
            table=[[a,b],[c,d]]
            if fisher[i]<(threshold):
                temp={"text": i,"size":math.log(1/fisher[i])}
                word_fisher.append(temp)
        f3.write(str(word_fisher)+"\n")
        
        matrix = [[0 for x in range(count)] for x in range(count)]
        f7.write(">"+ii+"\n")
        for i in clusterid:
            clusterid[i]=list(set(clusterid[i]))
            f4.write(i+"\t"+" ".join(clusterid[i])+"\n")
            f7.write(i+"\t"+"\t".join(clusterid[i])+"\n")
            for j in clusterid[i]:
                matrix[nmap[j]["id"]][nmap[i]["id"]]=matrix[nmap[j]["id"]][nmap[i]["id"]]+1
        for i in clusterdb:
            row=[0]*count
            f5.write(i+"\t"+" ".join(clusterdb[i])+"\n")
            for j in clusterdb[i]:
                row[nmap[j]["id"]]=nmap[i]["id"]
                matrix[nmap[j]["id"]][nmap[i]["id"]]=matrix[nmap[j]["id"]][nmap[i]["id"]]+1
        mat[ii]=(matrix,nmap)
    for i in mat:
        f2=open(folder+"/static/results/"+uuid+"/"+ids+"_graph/"+i+"_mat.txt","w")
        f2.write(str(mat[i][0])+"\n")
        f2.write(str(mat[i][1])+"\n")
        f9.write(i+"\t"+str(len(mat[i][1]))+"\n")
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    f6.close()
    f7.close()
    f8.close()
    f9.close()
    return col,nmap,domain_db,count_annotation,root_second_level,controls
