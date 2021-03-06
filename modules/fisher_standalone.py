#!/usr/bin/env python
# -*- coding: utf-8 -*-
#standalone
from gluon import *
import json
from scipy.stats import fisher_exact
import networkx as nx
import sys,os
def load(uuid,ids,protein,threshold,comp_list,folder,res_folder,path):
    fisher_annotation={}
    for ii in comp_list:
        fisher_annotation[ii]={}
        clusterdb={}
        clusterid={}
        prot_fisher=[]
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
        proteins_annotated={}
        while(seq!=""):
            seq=seq.strip().split("\t")
            fisher[seq[0]]=seq[1:]
            proteins_annotated[seq[0]]=""
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
        for i in protein:
            if fisher.has_key(i):
                for j in fisher[i]:
                    if fisherset.has_key(j):
                        annotation[j].append(i)
                        fisherset[j]=fisherset[j]+1
                    else:
                        annotation[j]=[]
                        annotation[j].append(i)
                        fisherset[j]=1
            else:
                fisherset_count.append(i)
        totalfisher=len(fisher)
        numberofproteins=len(protein)-len(list(set(fisherset_count)))
        fisher={}
        fisher_value={}
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
                fisher_annotation[ii][i]=""
                if fisher_value.has_key(fisher[i]):
                    fisher_value[fisher[i]].append(i)
                else:
                    fisher_value[fisher[i]]=[]
                    fisher_value[fisher[i]].append(i)
        f2=open(res_folder+"/"+ii+".txt","w")
        f2.write("id\tp_value\tproteins_involved\tdescription\tproteins\n")
        for i in fisher_value:
            for j in fisher_value[i]:
                f2.write(j+"\t"+str(i)+"\t"+str(len(annotation[j]))+"\t"+descr[j]+"\t"+" ".join(annotation[j])+"\n")
                fisher_json.append({"id":j,"p value":str(i),"proteins involved":str(len(annotation[j])),"description":descr[j],"proteins":annotation[j],"common":2})
        json.dump(fisher_json,open(folder+"/static/results/"+uuid+"/"+"/"+ids+"_fisher/"+ii+"fisher.json","w"))
    return fisher_annotation
