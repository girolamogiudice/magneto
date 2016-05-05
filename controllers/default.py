# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations
#default
#########################################################################
## This is a sample controller
## - index is the default action of any application
## - user is required for authentication and authorization
## - download is for downloading files uploaded in the db (does streaming)
#########################################################################

def index():
    return dict()

def test():
    return dict()
def test_table():
    return dict()

def tsfa():
    import os
    import types
    from types import NoneType
    from types import StringType
    from types import ListType
    expr_db={0:"protein_atlas_basal",1:"protein_atlas_cancer",2:"human_proteome_map_adult",3:"human_proteome_map_fetal",4:"PaxDb_integrated",5:"PeptideAatlas",6:"ProteomicsDB"}
    net_db={0:"intact",1:"intact_high_confidence",2:"biogrid"}
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
    PaxDb_integrated={0: 'brain', 2: 'colon', 3: 'esophagus', 4: 'female_gonad', 5: 'gall_bladder', 6: 'heart', 7: 'kidney', 8: 'liver', 9: 'lung', 10: 'pancreas', 11: 'placenta',
                12: 'plasma', 13: 'platelet', 14: 'prostate_gland', 15: 'rectum', 16: 'saliva', 17: 'skin', 18: 'testis', 19: 'urine', 20: 'uterus', 21: 'whole_organism'}

    db.define_table('input_net',
        Field('upload','upload',uploadfolder=os.path.join(request.folder,'upload'),requires=IS_NOT_EMPTY()),
        Field("graph_db",requires=IS_EMPTY_OR(IS_IN_SET(net_db,multiple=False))),
        Field('expr_db',requires=IS_EMPTY_OR(IS_IN_SET(expr_db,multiple=False))),
        Field('pa_basal',requires=IS_EMPTY_OR(IS_IN_SET(protein_atlas_basal,multiple=True))),
        Field('pa_cancer',requires=IS_EMPTY_OR(IS_IN_SET(protein_atlas_cancer,multiple=True))),
        Field('human_proteome_map_adult',requires=IS_EMPTY_OR(IS_IN_SET(human_proteome_map_adult,multiple=True))),
        Field('human_proteome_map_fetal',requires=IS_EMPTY_OR(IS_IN_SET(human_proteome_map_fetal,multiple=True))),
        Field('PaxDb_integrated',requires=IS_EMPTY_OR(IS_IN_SET(PaxDb_integrated,multiple=True))),
        Field('proteomicsdb',requires=IS_EMPTY_OR(IS_IN_SET(proteomicsdb,multiple=True))),
        Field('peptideatlas',requires=IS_EMPTY_OR(IS_IN_SET(peptideatlas,multiple=True))),
        Field("background","string",requires=IS_IN_SET(["Proteome","Human_uniprot"],multiple=False)),
        Field('Coexpression',requires=IS_EMPTY_OR(IS_IN_SET(["yes","no"],multiple=False))),
        Field('Show_advanced_options','boolean'),
        Field("threshold","float",requires=IS_FLOAT_IN_RANGE(0.0,0.1)),
        Field('Component_weight',requires=IS_FLOAT_IN_RANGE(0.0,1)),
        Field('Function_weight',requires=IS_FLOAT_IN_RANGE(0.0,1)),
        Field('Process_weight',requires=IS_FLOAT_IN_RANGE(0.0,1)),
        Field('Reactome_weight',requires=IS_FLOAT_IN_RANGE(0.0,1)),
        Field('Kegg_weight',requires=IS_FLOAT_IN_RANGE(0.0,1)))
        
    db.input_net.Coexpression.default = "yes"
    db.input_net.background.default = "Proteome"
    db.input_net.threshold.default = 0.05
    db.input_net.Component_weight.default=1.0
    db.input_net.Function_weight.default=1.0
    db.input_net.Process_weight.default=1.0
    db.input_net.Reactome_weight.default=1.0
    db.input_net.Kegg_weight.default=1.0
    db.input_net.threshold.show_if = (db.input_net.Show_advanced_options==True)
    db.input_net.Component_weight.show_if = (db.input_net.Show_advanced_options==True)
    db.input_net.Function_weight.show_if = (db.input_net.Show_advanced_options==True)
    db.input_net.Process_weight.show_if = (db.input_net.Show_advanced_options==True)
    db.input_net.Reactome_weight.show_if = (db.input_net.Show_advanced_options==True)
    db.input_net.Kegg_weight.show_if = (db.input_net.Show_advanced_options==True)
    db.input_net.pa_basal.show_if = (db.input_net.expr_db==0)
    db.input_net.pa_basal.hide_if = (db.input_net.expr_db!=0)
    db.input_net.pa_cancer.show_if = (db.input_net.expr_db==1)
    db.input_net.pa_cancer.hide_if = (db.input_net.expr_db!=1)
    db.input_net.human_proteome_map_adult.show_if = (db.input_net.expr_db==2)
    db.input_net.human_proteome_map_adult.hide_if = (db.input_net.expr_db!=2)
    db.input_net.human_proteome_map_fetal.show_if = (db.input_net.expr_db==3)
    db.input_net.human_proteome_map_fetal.hide_if = (db.input_net.expr_db!=3)
    db.input_net.PaxDb_integrated.show_if = (db.input_net.expr_db==4)
    db.input_net.PaxDb_integrated.hide_if = (db.input_net.expr_db!=4)
    db.input_net.proteomicsdb.show_if = (db.input_net.expr_db==6)
    db.input_net.proteomicsdb.hide_if = (db.input_net.expr_db!=6)
    db.input_net.peptideatlas.show_if = (db.input_net.expr_db==5)
    db.input_net.peptideatlas.hide_if = (db.input_net.expr_db!=5)
    db.input_net.Coexpression.widget = lambda field,value: SQLFORM.widgets.radio.widget(field,value,_style="text-align:left")
    db.input_net.background.widget = lambda field,value: SQLFORM.widgets.radio.widget(field,value,_style="text-align:left")
    db.input_net.pa_basal.widget = SQLFORM.widgets.multiple.widget
    db.input_net.pa_cancer.widget = SQLFORM.widgets.multiple.widget
    db.input_net.human_proteome_map_adult.widget = SQLFORM.widgets.multiple.widget
    db.input_net.human_proteome_map_fetal.widget = SQLFORM.widgets.multiple.widget
    db.input_net.PaxDb_integrated.widget = SQLFORM.widgets.multiple.widget
    db.input_net.proteomicsdb.widget = SQLFORM.widgets.multiple.widget
    db.input_net.peptideatlas.widget = SQLFORM.widgets.multiple.widget
    
    form = SQLFORM(db.input_net,_class="form-horizontal col-md-offset-3 col-md-4",labels = {'expr_db':'Protein Expression DB','threshold':'Significance Threshold' },_id="form_tsfa")
    if form.process().accepted:
        bg={"Proteome":"proteome","Human_uniprot":"uniprot"}
        from datetime import timedelta as timed
        import uuid,time
        #import filter.py
        uid=str(uuid.uuid4())
        if form.vars.pa_basal==None and form.vars.pa_cancer==None and form.vars.human_proteome_map_adult==None and form.vars.human_proteome_map_fetal==None and form.vars.PaxDb_integrated==None and form.vars.proteomicsdb==None and form.vars.peptideatlas==None:
            form.errors.protein_atlas_basal = "Cannot be empty!"
            response.flash = "Select at least a tissue!"
        if not isinstance( request.vars.pa_basal,NoneType) and isinstance(request.vars.pa_basal, StringType):
            request.vars.pa_basal=[request.vars.pa_basal]
        if not isinstance( request.vars.pa_cancer,NoneType) and isinstance(request.vars.pa_cancer, StringType):
            request.vars.pa_cancer=[request.vars.pa_cancer]
        if not isinstance( request.vars.human_proteome_map_adult,NoneType) and isinstance(request.vars.human_proteome_map_adult, StringType):
            request.vars.human_proteome_map_adult=[request.vars.human_proteome_map_adult]
        if not isinstance( request.vars.human_proteome_map_fetal,NoneType) and isinstance(request.vars.human_proteome_map_fetal, StringType):
            request.vars.human_proteome_map_fetal=[request.vars.human_proteome_map_fetal]
        if not isinstance( request.vars.PaxDb_integrated,NoneType) and isinstance(request.vars.PaxDb_integrated, StringType):
            request.vars.PaxDb_integrated=[request.vars.PaxDb_integrated]
        if not isinstance( request.vars.peptideatlas,NoneType) and isinstance(request.vars.peptideatlas, StringType):
            request.vars.peptideatlas=[request.vars.peptideatlas]
        if not isinstance( request.vars.proteomicsdb,NoneType) and isinstance(request.vars.proteomicsdb, StringType):
            request.vars.proteomicsdb=[request.vars.proteomicsdb]
        tissue_files=[]
        if request.vars.expr_db=="0":
             tissue_selection=request.vars.pa_basal
             for i in tissue_selection:
                tissue_files.append(protein_atlas_basal[int(i)])
        if request.vars.expr_db=="1":
             tissue_selection=request.vars.pa_cancer
             for i in tissue_selection:
                tissue_files.append(protein_atlas_cancer[int(i)])
        if request.vars.expr_db=="2":
             tissue_selection=request.vars.human_proteome_map_adult
             for i in tissue_selection:
                tissue_files.append(human_proteome_map_adult[int(i)])
        if request.vars.expr_db=="3":
             tissue_selection=request.vars.human_proteome_map_fetal
             for i in tissue_selection:
                tissue_files.append(human_proteome_map_fetal[int(i)])
        if request.vars.expr_db=="4":
             tissue_selection=request.vars.PaxDb_integrated
             for i in tissue_selection:
                tissue_files.append(PaxDb_integrated[int(i)])
        if request.vars.expr_db=="5":
             tissue_selection=request.vars.peptideatlas
             for i in tissue_selection:
                tissue_files.append(peptideatlas[int(i)])
        if request.vars.expr_db=="6":
             tissue_selection=request.vars.proteomicsdb
             for i in tissue_selection:
                tissue_files.append(proteomicsdb[int(i)])
        values=[uid,request.folder+"/upload/"+form.vars.upload,request.vars.graph_db,float(request.vars.threshold),request.vars.expr_db,
                tissue_files,form.vars.Coexpression,request.folder,float(request.vars.Component_weight),float(request.vars.Function_weight),float(request.vars.Process_weight),
                float(request.vars.Reactome_weight),float(request.vars.Kegg_weight),bg[form.vars.background]]
        
        if not os.path.exists(request.folder+"/static/results/"+uid):
            os.makedirs(request.folder+"/static/results/"+uid)
      
        
        db.query_tfsa.insert(uuid=uid,upload=request.vars.upload.value,upload_filename=form.vars.upload,threshold=form.vars.threshold,
                             graph_db=form.vars.graph_db,expr_db=form.vars.expr_db,protein_atlas_basal=form.vars.pa_basal,protein_atlas_cancer=form.vars.pa_cancer,
                             humanproteomemap_adult=form.vars.human_proteome_map_adult,humanproteomemap_fetal=form.vars.human_proteome_map_fetal,PaxDb_integrated=form.vars.PaxDb_integrated,
                             proteomicsdb=form.vars.proteomicsdb,peptideatlas=form.vars.peptideatlas,coexpression=form.vars.Coexpression,background=bg[form.vars.background],
                             Component_weight=form.vars.Component_weight,Function_weight=form.vars.Function_weight,Process_weight=form.vars.Process_weight,Reactome_weight=form.vars.Reactome_weight,
                             Kegg_weight=form.vars.Kegg_weight)
        #import mcn
        #mcn.tsfa(uid,request.folder+"/upload/"+form.vars.upload,request.vars.graph_db,float(request.vars.threshold),request.vars.expr_db,
        #       tissue_files,form.vars.Coexpression,request.folder,float(request.vars.Component_weight),float(request.vars.Function_weight),float(request.vars.Process_weight),
        #      float(request.vars.Reactome_weight),float(request.vars.Kegg_weight),bg[form.vars.background])
        s=scheduler.queue_task(tsfa,pargs=values,start_time=request.now,stop_time=request.now+timed(seconds=7200),timeout=3600)
        
        redirect(URL('wait_tsfa',args=[s.uuid,uid,request.vars.threshold,bg[form.vars.background],form.vars.graph_db]))
        
    return dict(form = form)
def wait_tsfa():
    import networkx as nx
    import fishertest,fisher_standalone_hierarchy
    
    import json
    net_db={"0":"intact","1":"intact_high_confidence","2":"biogrid"}
    sched_var=str(scheduler.task_status(request.args[0]).status)
    if scheduler.task_status(request.args[0]).status=="QUEUED" or scheduler.task_status(request.args[0]).status=="ASSIGNED":
        return dict(sched_var=sched_var)
    if scheduler.task_status(request.args[0]).status=="RUNNING":
        return dict(sched_var=sched_var)
    if scheduler.task_status(request.args[0]).status=="TIMEOUT" or scheduler.task_status(request.args[0]).status=="EXPIRED":
        redirect(URL('timeout'))
    if scheduler.task_status(request.args[0]).status=="FAILED":
        redirect(URL('failed'))

    if scheduler.task_status(request.args[0]).status=="COMPLETED":
        f1=open(request.folder+"static/results/"+request.args[1]+"/mapping_column.txt","r")
        seq=f1.readline()
        mapping={}
        while(seq!=""):
            seq=seq.strip().split("\t")
            mapping[seq[0]]=seq[1]
            seq=f1.readline()
        f1.close()
        
        choice=request.args[3]
        sample_number=len(mapping)
        
        descr_gene={}
        if choice=="proteome":
            path=request.folder+"/data/annotation/"+choice+"/"
        if choice=="uniprot":
            path=request.folder+"/data/annotation/"+choice+"/"
        #mapping protein_gene_name
        f1=open(path+"descr.txt")
        proteome_descr={}
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            proteome_descr[seq[0]]=seq[1]
            seq=f1.readline()
        f1.close()
        
        protein_descr={}
        root={}
        root = {"name": "database","children": []}
        root_second_level={}
        count_annotation={}
        count_protein={}
        db_hierarchy={'DB':['Approved', 'Small_molecule', 'Experimental', 'Nutraceutical', 'Illicit', 'Withdrawn', 'Investigational', 'Biotech'],
                      'KDi':['Respiratory_diseases', 'Digestive_system_diseases', 'Congenital_disorders_of_metabolism', 'Urinary_system_diseases', 'Musculoskeletal_diseases', 'Cardiovascular_diseases',
                             'Immune_system_diseases', 'Cancers', 'Endocrine_and_metabolic_diseases', 'Reproductive_system_diseases', 'Skin_diseases', 'Nervous_system_diseases', 'Other_congenital_disorders'],
                      'KDr':['Hormonal_Agents_Suppressant_(Parathyroid)', 'Antiparkinson_Agents', 'Therapeutic_Nutrients_Minerals_Electrolytes', 'Sleep_Disorder_Agents', 'Antifungals', 'Analgesics',
                             'Hormonal_Agents_Stimulant_Replacement_Modifying_(Adrenal)', 'Central_Nervous_System_Agents', 'Anti_inflammatory_Agents', 'Bipolar_Agents', 'Metabolic_Bone_Disease_Agents',
                             'Dental_and_Oral_Agents', 'Antidementia_Agents', 'Gastrointestinal_Agents', 'Antimyasthenic_Agents', 'Antivirals', 'Otic_Agents', 'Antispasticity_Agents',
                             'Hormonal_Agents_Stimulant_Replacement_Modifying_(Pituitary)', 'Genitourinary_Agents', 'Blood_Products_Modifiers_Volume_Expanders', 'Antigout_Agents', 
                             'Respiratory_Tract_Pulmonary_Agents', 'Hormonal_Agents_Stimulant_Replacement_Modifying_(Prostaglandins)', 'Hormonal_Agents_Suppressant_(Thyroid)',
                             'Hormonal_Agents_Stimulant_Replacement_Modifying_(Thyroid)', 'Hormonal_Agents_Stimulant_Replacement_Modifying_(Sex_Hormones_Modifiers)', 'Antiemetics', 'Cardiovascular_Agents',
                             'Anesthetics', 'Anticonvulsants', 'Immunological_Agents', 'Antibacterials', 'Inflammatory_Bowel_Disease_Agents', 'Dermatological_Agents',
                             'Anti_Addiction_Substance_Abuse_Treatment_Agents', 'Hormonal_Agents_Suppressant_(Pituitary)', 'Hormonal_Agents_Suppressant_(Adrenal)',
                             'Antimycobacterials', 'Antipsychotics', 'Antineoplastics', 'Enzyme_Replacement_Modifiers', 'Blood_Glucose_Regulators', 'Antiparasitics',
                             'Anxiolytics', 'Antimigraine_Agents', 'Skeletal_Muscle_Relaxants', 'Ophthalmic_Agents', 'Antidepressants'],
                      'K':['Aging', 'Target_based_classification_G_protein_coupled_receptors', 'Cancers_Specific_types', 'Endocrine_system', 'Transcription', 'Energy_metabolism', 'Cellular_community',
                           'Substance_dependence', 'Sensory_system', 'Chemical_structure_transformation_maps', 'Structure_based_classification', 'Nervous_system', 'Circulatory_system',
                           'Signaling_molecules_and_interaction', 'Nucleotide_metabolism', 'Metabolism_of_terpenoids_and_polyketides', 'Chronology_Other_drugs', 'Skeleton_based_classification',
                           'Infectious_diseases_Bacterial', 'Amino_acid_metabolism', 'Immune_diseases', 'Transport_and_catabolism', 'Translation', 'Chronology_Nervous_system_agents', 'Cell_motility',
                           'Replication_and_repair', 'Endocrine_and_metabolic_diseases', 'Development', 'Drug_resistance','Digestive_system', 'Cell_growth_and_death', 'Signal_transduction', 'Cancers_Overview',
                           'Cardiovascular_diseases', 'Xenobiotics_biodegradation_and_metabolism', 'Environmental_adaptation', 'Chronology_Antiinfectives', 'Target_based_classification_Transporters',
                           'Carbohydrate_metabolism', 'Biosynthesis_of_other_secondary_metabolites', 'Target_based_classification_Enzymes', 'Glycan_biosynthesis_and_metabolism', 'Lipid_metabolism',
                           'Excretory_system', 'Metabolism_of_other_amino_acids', 'Metabolism_of_cofactors_and_vitamins', 'Target_based_classification_Ion_channels', 'Membrane_transport',
                           'Target_based_classification_Nuclear_receptors', 'Infectious_diseases_Viral', 'Neurodegenerative_diseases', 'Folding_sorting_and_degradation', 'Immune_system',
                           'Global_and_overview_maps', 'Chronology_Antineoplastics', 'Infectious_diseases_Parasitic'],
                      'Or':['Rare_neurological_diseases', 'Rare_abdominal_surgical_diseases', 'Rare_odontological_diseases', 'Rare_systemic_and_rhumatological_diseases', 'Rare_urogenital_diseases',
                            'Rare_cardiac_diseases', 'Rare_genetic_diseases', 'Rare_hepatic_diseases', 'Rare_intoxications', 'Rare_respiratory_diseases', 'Developmental_anomalies_during_embryogenesis',
                            'Rare_renal_diseases', 'Rare_immunological_diseases', 'Rare_haematological_diseases', 'Rare_gastroenterological_diseases', 'Rare_tumors', 'Rare_allergic_disease',
                            'Inborn_errors_of_metabolism', 'Teratologic_disorders', 'Rare_infectious_diseases', 'Rare_skin_diseases', 'Rare_endocrine_diseases', 'Rare_infertility_disorders',
                            'Rare_cardiac_malformations', 'Rare_surgical_maxillo-facial_diseases', 'Rare_otorhinolaryngological_diseases', 'Rare_bone_diseases', 'Rare_eye_diseases',
                            'Rare_circulatory_system_diseases', 'Rare_gynaecological_and_obstetric_diseases', 'Rare_surgical_thoracic_diseases'],
                      'R':['Cell_Cell_communication', 'Cell_Cycle', 'Cellular_responses_to_stress', 'Chromatin_organization', 'Circadian_Clock', 'DNA_Repair', 'DNA_Replication', 'Developmental_Biology',
                           'Disease', 'Extracellular_matrix_organization', 'Gene_Expression', 'Hemostasis', 'Immune_System', 'Metabolism', 'Metabolism_of_proteins', 'Mitophagy', 'Muscle_contraction',
                           'Neuronal_System', 'Organelle_biogenesis_and_maintenance', 'Programmed_Cell_Death', 'Reproduction', 'Signal_Transduction', 'Transmembrane_transport_of_small_molecules',
                           'Vesicle_mediated_transport']}


        fisher_annotation={}
        fisher_annotation_hierarchy={}
        domain_db={}
        temp_db={}
        children={}
        alias_dict={"C":"Component","F":"Function","P":"Process","K":"Kegg Pathway","R":"Reactome","O":"Omim","KDi":"Kegg Disease","KDr":"Kegg Drug","DB":"Drugbank","Or":"Orphanet",
                    "HPi":"Virus","T":"toxins"}
        databases={"Gene Ontology":["C","F","P"],"Pathway":["R","K"],"Disease":["O","Or","KDi"],"Drugs":["DB","KDr"],"Toxins":["T"],"Virus":["HPi"]}
        for i in ["C","P","F","R","K","O","KDr","KDi","DB","Or","HPi","T"]:
            fisher_annotation[i]=[]
            domain_db[i]=[]
            count_annotation[i]={}
        for i in ["R","DB","K","KDi","KDr","Or"]:
            fisher_annotation_hierarchy[i]={}
            
        nodes_in_the_graph={}
        f1=open(request.folder+"static/results/"+request.args[1]+"/nodes_present_in_graph.txt","r")
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            nodes_in_the_graph[seq[0]]=seq[1:]
            seq=f1.readline()

        f1=open(request.folder+"static/results/"+request.args[1]+"/nodes.txt","r")
        seq=f1.readline()
        fisher_alone={}
        fisher_alone_hierarchy={}
        start_nodes={}
        protein_stats={}
        temp=[]
        not_annotated={}
        proteome_annotated={}
        db_background=request.args[3] 
        while(seq!=""):
            seq=seq.strip().split("\t")
            fisher_alone[seq[0]]={}
            fisher_alone_hierarchy[seq[0]]={}
            start_nodes[seq[0]]=seq[1:]
            not_annotated[seq[0]]=[]
            proteome_annotated[seq[0]]=[]
            protein_stats[seq[0]]={}
            fisher_alone[seq[0]]=fisher_standalone.load(request.args[1],seq[1:],float(request.args[2]),["C","P","F","R","K","O","KDr","KDi","DB","Or","T","HPi"],request.folder,
                               request.folder+"/static/results/"+request.args[1]+"/"+seq[0]+"_fisher",path,fisher_annotation)
            fisher_alone_hierarchy[seq[0]]=fisher_standalone_hierarchy.load(request.args[1],seq[1:],float(request.args[2]),db_hierarchy,request.folder,
                               request.folder+"/static/results/"+request.args[1]+"/"+seq[0]+"_fisher",path,fisher_annotation_hierarchy,db_background)
            
            seq=f1.readline()
        f1=open(request.folder+"static/results/"+request.args[1]+"/nodes_graph.txt","r")
        nodes_graph={}
        #proteome whole proteome e le network sono diverse
        seq=f1.readline()
        
        while(seq!=""):
            sample=mapping[seq[0]]
            f2=open(request.folder+"static/results/"+request.args[1]+"/"+seq[0]+"_graph/not_present.txt","r")
            not_present_in_network=f2.readline().strip().split("\t")
            f2.close()
            seq=seq.strip().split("\t")
            
            for j in seq[1:]:
                if proteome_descr.has_key(j):
                    protein_descr[j]=proteome_descr[j]
                    
                else:
                    protein_descr[j]="N/A"
                if protein_stats[seq[0]].has_key(j):
                    protein_stats[seq[0]][j].append(sample)
                else:
                    protein_stats[seq[0]][j]=[]
                    protein_stats[seq[0]][j].append(sample)
            for j in not_present_in_network:
                if proteome_descr.has_key(j):
                    proteome_annotated[seq[0]].append(j)
                    protein_descr[j]=proteome_descr[j]
                else:
                    not_annotated[seq[0]].append(j)
                if protein_stats[seq[0]].has_key(j):
                    protein_stats[seq[0]][j].append(sample)
                else:
                    protein_stats[seq[0]][j]=[]
                    protein_stats[seq[0]][j].append(sample)
            #devo controllare il fisher test e sbagliato check the path        
            col,nmap,domain_db,count_annotation,root_second_level=fishertest.load(request.args[1],seq[0],seq[1:]+proteome_annotated[seq[0]],
                                                                                       float(request.args[2]),["C","P","F","R","K","O","KDr","KDi","DB","Or","T","HPi"],
                            request.args[3],request.folder,start_nodes[seq[0]],domain_db,protein_descr,count_annotation,root_second_level,fisher_alone[seq[0]],sample,sample_number)
            
            fisher_alone_hierarchy[seq[0]]=fisher_standalone_hierarchy.load(request.args[1],seq[1:]+proteome_annotated[seq[0]],float(request.args[2]),db_hierarchy,request.folder,
                               request.folder+"/static/results/"+request.args[1]+"/"+seq[0]+"_graph",path,fisher_annotation_hierarchy,db_background)
                        
            nodes_graph[seq[0]]=seq[1:]
            graph_mcn=nx.read_gpickle(request.folder+"static/results/"+request.args[1]+"/"+seq[0]+"_graph/graph_mcn.gpickle")
            number_of_nodes=graph_mcn.number_of_nodes()
            count=0
            nodes=[]
            edges=[]
            
            f2=open(request.folder+"static/results/"+request.args[1]+"/"+seq[0]+"_graph/prot_descr.txt","w")
            for i in graph_mcn.nodes():
                graph_mcn.node[i]['id']=i
                f2.write(i+"\t"+protein_descr[i]+"\n")
                graph_mcn.node[i]['label']=i+"("+protein_descr[i]+")"
                graph_mcn.node[i]['size']=4
                graph_mcn.node[i]['x']=math.cos(2 * count * 3.14 / number_of_nodes)
                graph_mcn.node[i]['y']=math.sin(2 * count * 3.14 / number_of_nodes)
                graph_mcn.node[i]['color']=['#6495ED']
                if col.has_key(i):
                    graph_mcn.node[i]['colors']=list(set(col[i]))
                else:
                    graph_mcn.node[i]['colors']=['#ABFFE3']

                nodes.append(graph_mcn.node[i])
                count=count+1
            for i in nodes_in_the_graph[seq[0]]:
                graph_mcn.node[i]['color']=['#FF0000']
            count=0
            f2.close()
            for i in graph_mcn.edges():
                graph_mcn.edge[i[0]][i[1]]['id']="e"+str(count)
                graph_mcn.edge[i[0]][i[1]]['target']=i[1]
                graph_mcn.edge[i[0]][i[1]]['source']=i[0]
                graph_mcn.edge[i[0]][i[1]]['size']=7
                graph_mcn.edge[i[0]][i[1]]['color']='#C0C0C0'
                count=count+1
                edges.append(graph_mcn.edge[i[0]][i[1]])
            json.dump(dict(nodes=nodes,edges=edges),open(request.folder+"static/results/"+request.args[1]+"/"+seq[0]+"_graph/graph_mcn.json",'w'))
            seq=f1.readline()

        if sample_number>1:
            for i in databases:
                flag=0
                for j in databases[i]:
                    if len(count_annotation[j].values())>0:
                        if flag==0:
                            root["children"].append({"name":i,"children":[]})
                            flag=1
                        if len(count_annotation[j])!=0:
                            root["children"][-1]["children"].append({"name":alias_dict[j],"children":count_annotation[j].values()})

            json.dump(root,open(request.folder+"static/results/"+request.args[1]+"/stats.json",'w'))


        for i in ["C","P","F","R","K","O","KDr","KDi","DB","Or","T","HPi"]:
            f2=open(request.folder+"/static/results/"+request.args[1]+"/"+i+"_domain.txt","w")
            f2.write(str(list(set(domain_db[i]))))
            f2.close()
        return dict(sched_var=sched_var)
       
    return dict(sched_var=sched_var)

def results_tsfa():
    f1=open(request.folder+"static/results/"+request.args[0]+"/mapping_column.txt","r")
    seq=f1.readline()
    mapping_column=[]
    mapping={}
    while(seq!=""):
        seq=seq.strip().split("\t")
        mapping_column.append((seq[0],seq[1].strip()))
        mapping[seq[1].strip()]=seq[0]
        seq=f1.readline()
    sample_number=len(mapping)
    f1.close()
    f1=open(request.folder+"static/results/"+request.args[0]+"/nodes.txt","r")
    starting_nodes={}
    seq=f1.readline()
    while(seq!=""):
        seq=seq.strip().split("\t")
        starting_nodes[seq[0]]=seq[1:]
        seq=f1.readline()
    f1.close()
    menu_list={}

    alias_dict={"C":"Component","F":"Function","P":"Process","K":"Kegg Pathway","R":"Reactome","O":"Omim","KDi":"Kegg Disease","KDr":"Kegg Drug","DB":"Drugbank","Or":"Orphanet",
                "T":"Toxins","HPi":"Virus"}

    menu_mapping={}
    for i in mapping_column:
        menu_list[i[0]]={}
        menu_mapping[i[0]]={}
        for j in ["R","K","KDr","KDi","Or","DB"]:
            menu_mapping[i[0]][j]={}
            f1=open(request.folder+"static/results/"+request.args[0]+"/"+i[0]+"_graph/"+j+"_hierarchy/files_domain.txt","r")
            seq=f1.readline()
            menu_list[i[0]][j]=""
            #<ul class="dropdown-menu dm" id="Reactome" style="position:absolute;display: none;left:100%;top:0px;margin-top: 0px;margin-left:-5px;margin:-1px; padding:0;">  </ul>
            seq=seq.strip().split("\t")
            if len(seq)==1:
                continue
            count=0
            for k in seq:
                menu_list[i[0]][j]=  menu_list[i[0]][j]+ '<li class="test" id="'+k+'"><a tabindex="-1" href="#">'+k.replace("_"," ")+'</a><li>'
                menu_mapping[i[0]][j][k]=str(count)
                count=count+1
    databases={"Gene Ontology":["C","F","P"],"Pathway":["R","K"],"Disease":["O","Or","KDi"],"Drugs":["DB","KDr"],"Toxins":["T"],"Virus":["HPi"]}

    alias_colors={"C":"#d5f0f4","P":"#2171b5","F":"#6baed6","R":"#abe16c","K":"#5f8726","O":"#7d6396","KDr":"#b23131","KDi":"#ffcdff","DB":"#ff0000","Or":"#c19bce","HPi":"#A9A9A9","T":"#b56b19"}
    return dict(mapping_column=mapping_column,mapping=mapping,starting_nodes=starting_nodes,databases=databases,alias_colors=alias_colors,alias_dict=alias_dict,sample_number=sample_number,menu_list=menu_list,
               menu_mapping=menu_mapping)

def load_request():
    alias_dict={"C":"Component","F":"Function","P":"Process","K":"Kegg Pathway","R":"Reactome","O":"Omim","KDi":"Kegg Disease","KDr":"Kegg Drug","DB":"Drugbank","Or":"Orphanet","T":"Toxins","HPi":"Virus"}
    directory=request.args[0]
    ids=request.args[1]
    f1=open(request.folder+"static/results/"+request.args[0]+"/"+ids+"graph/not_present.txt")
    not_present=f1.readline()
    clusterid={}
    f1=open(request.folder+"static/results/"+request.args[0]+"/"+ids+"graph/graph_mcn.json","rU")
    data=f1.readline()
    f1.close()
    f1=open(request.folder+"static/results/"+request.args[0]+"/"+ids+"graph/cluster_id.txt")
    seq=f1.readline()
    while(seq!=""):
        if seq[0]==">":
            key=seq[1:].strip()
            clusterid[key]={}
        else:
            seq=seq.strip().split("\t")
            clusterid[key][seq[0]]=seq[1:]
        seq=f1.readline()
    neigh={}
    f1=open(request.folder+"/static/results/"+request.args[0]+"/"+ids+"graph/neigh.txt",'rU')
    seq=f1.readline()
    while(seq!=""):
        seq=seq.strip().split("\t")
        neigh[seq[0]]=seq[1:]
        seq=f1.readline()
    descr={}
    for i in ["C","P","F","R","K","O","KDr","KDi","DB","Or","T","HPi"]:
        descr[i]={}
        f1=open(request.folder+"/static/results/"+request.args[0]+"/"+ids+"graph/"+i+"_descr.txt",'rU')
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            descr[i][seq[0]]=seq[1]
            seq=f1.readline()
    protein_descr={}
    f2=open(request.folder+"static/results/"+request.args[0]+"/"+ids+"graph/prot_descr.txt","rU")
    seq=f2.readline()
    while(seq!=""):
        seq=seq.strip().split("\t")
       
        protein_descr[seq[0]]=seq[1]
        seq=f2.readline()
    return dict(directory=directory,data=data,clusterid=clusterid,neigh=neigh,descr=descr,protein_descr=protein_descr,alias_dict=alias_dict,not_present=not_present)

def load_stats():
    sample_number= request.args[2].strip("_")
    ids=request.args[1]
    f1=open(request.folder+"static/results/"+request.args[0]+"/"+ids+"_graph/stats.txt","r")
    stats=[]
    seq=f1.readline()
    while(seq!=""):
        stats.append(seq.strip())
        seq=f1.readline()

    path="/magneto/static/results/"+request.args[0]+"/stats.json"
    
    return dict(stats=stats,path=path,sample_number=sample_number)

def load_chord():
    values=request.env.http_web2py_component_location.split("#")
    ids=values[1]
    db_to_load=values[2]
    f1=open(request.folder+"static/results/"+request.args[0][0:-1]+"/"+ids+"_graph/"+db_to_load+"_mat.txt")
    mat=f1.readline()
    mmap=f1.readline()
    f1=open(request.folder+"static/results/"+request.args[0][0:-1]+"/"+ids+"_graph/"+db_to_load+"_cluster_id.txt","r")
    cluster={}
    seq=f1.readline()
    while(seq!=""):
        seq=seq.strip().split("\t")
        cluster[seq[0]]=seq[1]
        seq=f1.readline()
    f1=open(request.folder+"static/results/"+request.args[0][0:-1]+"/"+ids+"_graph/"+db_to_load+"_cluster_db.txt","r")
    seq=f1.readline()
    while(seq!=""):
        seq=seq.strip().split("\t")
        cluster[seq[0]]=seq[1]
        seq=f1.readline()
    return dict(mmap=mmap,mat=mat,db_to_load=db_to_load,cluster=cluster)

def load_treemap():
    values=request.env.http_web2py_component_location.split("#")
    ids=values[1]+"_"
    db_to_load=values[2]
    path="/magneto/static/results/"+request.args[0][0:-1]+"/"+ids+"graph/"+db_to_load+"fisher.txt"
    return dict(path=path)

def load_bubble():
    values=request.env.http_web2py_component_location.split("#")
    ids=values[1]
    db_to_load=values[2]
    path="/magneto/static/results/"+request.args[0][0:-1]+"/"+ids+"_graph/"+db_to_load+"fisher.txt"
    return dict(path=path)

def load_word_cloud():
    values=request.env.http_web2py_component_location.split("#")
    ids=values[1]
    db_to_load=values[2]
    path=request.folder+"static/results/"+request.args[0][0:-1]+"/"+ids+"_graph/"+db_to_load+"_word_fisher.txt"
    f1=open(path,"r")
    word_cloud={}
    word_cloud=f1.readline()
    return dict(word_cloud=word_cloud)
    
def tfa():
    import types,os
    from types import NoneType
    from types import StringType
    from types import ListType
    score={0:"bwmod",1:"degree",2:"bw",3:"pagerank"}
    net_db={0:"intact",1:"intact_hc",2:"biogrid"}
    db.define_table('input_net',
        Field('upload','upload',uploadfolder=os.path.join(request.folder,'upload'),requires=IS_NOT_EMPTY()),
        Field("graph_db",requires=IS_EMPTY_OR(IS_IN_SET(net_db,multiple=False))),
        Field("threshold","float",requires=IS_FLOAT_IN_RANGE(0.0,0.1)),
        Field('score',requires=IS_EMPTY_OR(IS_IN_SET(score,multiple=False))),
        Field("topology_ranking","text",requires=IS_NOT_EMPTY()))
    
    db.input_net.threshold.default = 0.05
    form = SQLFORM(db.input_net,_class="col-md-4 col-md-offset-2",_id="form_tfa",col3 = {'score':XML('<input type="button" class="btn btn-primary" value="Add filter" onClick="addtext();"><input type="button" class="btn btn-warning" value="Remove filters" onClick="removetext();">')}).process()
    submit = form.element("input",_type="submit")
    submit["_onclick"] = 'document.getElementById("input_net_topology_ranking").disabled = false;'
    if form.accepted:
        from datetime import timedelta as timed
        import uuid,time
        #import filter.py
        uid=str(uuid.uuid4())
        values=[uid,form.vars.upload,request.vars.graph_db,float(request.vars.threshold),form.vars.topology_ranking.split()]
        db.query_tfa.insert(uuid=uid,upload=request.vars.upload.value,upload_filename=form.vars.upload,threshold=form.vars.threshold,
                            graph_db=form.vars.graph_db,topology_ranking=form.vars.topology_ranking.split())
        s=scheduler.queue_task(tfa,pargs=values,start_time=request.now,stop_time=request.now+timed(seconds=7200),timeout=240)
        redirect(URL('wait_tfa',args=[s.uuid,uid]))
    return dict(form=form)

def wait_tfa():
    sched_var=str(scheduler.task_status(request.args[0]).status)
    if scheduler.task_status(request.args[0]).status=="QUEUED" or scheduler.task_status(request.args[0]).status=="ASSIGNED":
        return dict(sched_var=sched_var)
    if scheduler.task_status(request.args[0]).status=="RUNNING":
        return dict(sched_var=sched_var)
    if scheduler.task_status(request.args[0]).status=="TIMEOUT" or scheduler.task_status(request.args[0]).status=="EXPIRED":
        redirect(URL('timeout'))
    if scheduler.task_status(request.args[0]).status=="FAILED":
        redirect(URL('failed'))

    if scheduler.task_status(request.args[0]).status=="COMPLETED":
        return dict(sched_var=sched_var)

    return dict(sched_var=sched_var)

def results_tfa():
    return dict()

def help():
    return dict()

def retrieve():
    db.define_table('retrieve',
        Field("ids",requires = IS_ALPHANUMERIC(error_message='must be alphanumeric!')))
    form = SQLFORM(db.retrieve,_class="col-md-4 col-md-offset-2",_id="form_tfa").process()
    if form.accepted:
        if len(form.vars.ids)!=16:
            form.errors.ids = "The id should be 16 characters long, please check!"
            response.flash = "error!"
        redirect(URL('results_tfa',args=form.vars.ids))
    return dict(form=form)

def team():
    return dict()

def statistics():
    return dict()

def user():
    """
    exposes:
    http://..../[app]/default/user/login
    http://..../[app]/default/user/logout
    http://..../[app]/default/user/register
    http://..../[app]/default/user/profile
    http://..../[app]/default/user/retrieve_password
    http://..../[app]/default/user/change_password
    http://..../[app]/default/user/manage_users (requires membership in
    http://..../[app]/default/user/bulk_register
    use @auth.requires_login()
        @auth.requires_membership('group name')
        @auth.requires_permission('read','table name',record_id)
    to decorate functions that need access control
    """
    return dict(form=auth())


@cache.action()
def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request, db)


def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    return service()
