#!/usr/bin/env python
# -*- coding: utf-8 -*-
from gluon import *
from scipy.stats import fisher_exact
import sys,os,json
def load(uuid,protein,threshold,comp_list,folder,res_folder,path,fisher_alone_hierarchy,db_background,control,controls_h):
    databases=[]
    alias_dict={"K":"Kegg Pathway","R":"Reactome","KDi":"Kegg Disease","KDr":"Kegg Drug","DB":"Drugbank","Or":"Orphanet"}
    fisher_hierarchy={}
    for ii in comp_list:
        json_tree={}
        fisher_hierarchy={}
        
        json_tree[ii]={"name":alias_dict[ii],"children":[]}
        if not os.path.exists(res_folder+"/"+ii+"_hierarchy"):
            os.makedirs(res_folder+"/"+ii+"_hierarchy")
        descr={}
        sub_categories={}
        sub_categories[ii]=""
        
        f1=open(path+ii+"_descr.txt","r")
        seq=f1.readline()
        while (seq!=""):
            seq=seq.strip().split("\t")
            descr[seq[0]]=seq[1]
            seq=f1.readline()
        f1.close()

        for jj in comp_list[ii]:
            temp_table=[]
            controls_h[ii][jj]={}
            fisher_hierarchy[jj]={}
            prot_fisher=[]
            f1=open(path+ii+"_hierarchy_"+db_background+"/"+jj+".txt","rU")
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
            f1=open(path+ii+"_hierarchy_"+db_background+"/"+jj+"_count.txt","r")
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
            fisher_json=[]
            for i in fisherset:
                a=fisherset[i]
                b=numberofproteins-a
                c=fisher_count[i]-a
                d=totalfisher-a-b-c
                table=[[a,b],[c,d]]
                fisher[i]=fisher_exact(table,alternative ="greater")[1]
                if fisher[i]<(threshold/lenfisherset):
                    flag=1
                    controls_h[ii][jj][i]=""
                    #fisher_annotation[ii][i]=""
                    if fisher_value.has_key(fisher[i]):
                        fisher_value[fisher[i]].append(i)
                    else:
                        fisher_value[fisher[i]]=[]
                        fisher_value[fisher[i]].append(i)
            if flag==1:
                f2=open(res_folder+"/"+ii+"_hierarchy/"+jj+".txt","w")
                f2.write("id\tp_value\tproteins_involved\tdescription\tproteins\n")
                databases.append(ii)
                sub_categories[ii]=sub_categories[ii]+jj+"\t"
                temp=[]
                temp_alone=[]
                for i in fisher_value:
                    if len(fisher_value)>0:
                        for j in fisher_value[i]:
                            if fisher_alone_hierarchy[ii].has_key(jj):
                                if fisher_alone_hierarchy[ii][jj].has_key(j):
                                    temp_table.append({"id":j,"description":descr[j],"proteins":annotation[j],"common":"1","p value":str(i),"proteins involved":str(len(annotation[j]))})
                                    temp.append({"name":j,"description":descr[j],"value":str(i),"common":"1","proteins involved":" ".join(annotation[j])})
                                    del fisher_alone_hierarchy[ii][jj][j]
                                    if len(fisher_alone_hierarchy[ii][jj])==0:
                                        del fisher_alone_hierarchy[ii][jj]
                                else:
                                    temp_table.append({"id":j,"description":descr[j],"proteins":annotation[j],"common":"0","p value":str(i),"proteins involved":str(len(annotation[j]))})
                                    temp.append({"name":j,"description":descr[j],"value":str(i),"common":"0","proteins involved":" ".join(annotation[j])})
                            else:
                                temp_table.append({"id":j,"description":descr[j],"proteins":annotation[j],"common":"0","p value":str(i),"proteins involved":str(len(annotation[j]))})
                                temp.append({"name":j,"description":descr[j],"value":str(i),"common":"0","proteins involved":" ".join(annotation[j])})
                            f2.write(j+"\t"+str(i)+"\t"+descr[j]+"\t"+" ".join(annotation[j])+"\n")
                
                if fisher_alone_hierarchy[ii].has_key(jj):
                    
                    for k in fisher_alone_hierarchy[ii][jj]:
                        temp.append({"name":k,"description":descr[k],"value":"","proteins involved":"","common":"2"})
                        
                json_tree[ii]["children"].append({"name":jj,"children":temp})
                f2.close()

            if fisher_alone_hierarchy[ii].has_key(jj):
                if flag==0:
                    sub_categories[ii]=sub_categories[ii]+jj+"*"+"\t"
                    temp=[]
                    for k in fisher_alone_hierarchy[ii][jj]:
                        temp.append({"name":k,"description":descr[k],"value":"","proteins involved":"","common":"2"})
                    json_tree[ii]["children"].append({"name":jj,"children":temp})

                for i in fisher_alone_hierarchy[ii][jj]:
                    temp_table.append({"id":i,"description":descr[i],"proteins":"","common":"2","p value":"-1","proteins involved":"0"})

                
            if len(temp_table)>=1:
                json.dump(temp_table,open(res_folder+"/"+ii+"_hierarchy/"+jj+"fisher.json","w"))
        
        if len(json_tree[ii]["children"])>0:
            json.dump(json_tree[ii],open(res_folder+"/"+ii+"_hierarchy/json.txt","w"))

        for i in fisher_hierarchy[jj]:
            if fisher_annotation[ii].has_key(i):
                fisher_json.append({"id":i,"p value":str(fisher_value[i]),"proteins involved":str(len(annotation[i])),"description":descr[i],"proteins":annotation_gene[i],"common":1})
            elif not fisher_value.has_key(i):
                fisher_json.append({"id":i,"p value":"-1","proteins involved":str(len(fisher_annotation[ii][i])),"description":descr[i],"common":2})
            else:
                fisher_json.append({"id":i,"p value":str(fisher_value[i]),"proteins involved":str(len(annotation[i])),"description":descr[i],"proteins":annotation_gene[i],"common":0})
        
        json.dump(fisher_json,open(res_folder+"/"+ii+"_hierarchy/"+jj+"fisher.json","w"))
        f2=open(res_folder+"/"+ii+"_hierarchy/files_domain.txt","w")
        f2.write(sub_categories[ii].strip())
        f2.close()
        f2=open(res_folder+"/hierarchy_databases.txt","w")
        f2.write("\t".join(list(set(databases))))
        f2.close()
    return fisher_hierarchy,controls_h
