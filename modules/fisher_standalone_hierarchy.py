#!/usr/bin/env python
# -*- coding: utf-8 -*-
from gluon import *
from scipy.stats import fisher_exact
import sys,os
def load(uuid,protein,threshold,comp_list,folder,res_folder,path,fisher_annotation):
    for ii in comp_list:
        if not os.path.exists(res_folder+"/"+ii+"_hierarchy"):
            os.makedirs(res_folder+"/"+ii+"_hierarchy")
        descr={}
        f1=open(path+ii+"_descr.txt","r")
        seq=f1.readline()
        while (seq!=""):
            seq=seq.strip().split("\t")
            descr[seq[0]]=seq[1]
            seq=f1.readline()
        f1.close()
        for jj in comp_list[ii]:
            
            prot_fisher=[]
            f1=open(path+ii+"_hierarchy/"+jj+".txt","r")
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
            f1=open(path+ii+"_hierarchy/"+jj+"_count.txt","r")
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
            flag=0
            for i in fisherset:
                a=fisherset[i]
                b=numberofproteins-a
                c=fisher_count[i]-a
                d=totalfisher-a-b-c
                table=[[a,b],[c,d]]
                fisher[i]=fisher_exact(table,alternative ="greater")[1]
                if fisher[i]<(threshold/lenfisherset):
                    flag=1
                    fisher_annotation[ii][i]=""
                    if fisher_value.has_key(fisher[i]):
                        fisher_value[fisher[i]].append(i)
                    else:
                        fisher_value[fisher[i]]=[]
                        fisher_value[fisher[i]].append(i)
            if flag==1:
                f2=open(res_folder+"/"+ii+"_hierarchy/"+jj+".txt","w")
                for i in fisher_value:
                    if len(fisher_value)>0:
                        for j in fisher_value[i]:
                            f2.write(j+"\t"+str(i)+"\t"+descr[j]+"\t"+" ".join(annotation[j])+"\n")
    return fisher_annotation
