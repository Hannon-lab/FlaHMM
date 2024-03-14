import os
import sklearn

import pandas as pd
import numpy as np
from hmmlearn import hmm
try: # version > 0.2.7
   from hmmlearn.hmm import CategoricalHMM as MultinomialHMM
except: # version <= 0.2.7
   from hmmlearn.hmm import MultinomialHMM
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns

#Read Coordinates excel file, remove NaN rows and keep the columns with coordinates for flam
#Change '+','-' signs to 'plus' and 'minus' respectively
coordinates=pd.read_excel('data/TableS1.xlsx')
coordinates=coordinates[coordinates['Region'].isin(['flam','flam-5p'])].reset_index(drop=True) # Prefer 5' end over 3' end
coordinates['Strand']=np.where(coordinates['Strand']=='+', 'plus', 'minus')

#Create a Dictionary, where key is the species with the build id and the values are the coordinates for the flam 
build2coords={}
species2build={}
for row in coordinates.iterrows():
    if row[1]['Species'] not in species2build.keys():
        species2build[row[1]['Species']]=[row[1]['Assembly']]
    else:
        update_build_list=species2build[row[1]['Species']]+[row[1]['Assembly']]
        species2build[row[1]['Species']]=update_build_list


for row in coordinates.iterrows():
    build2coords[row[1]['Species']+'_build_'+row[1]['Assembly']]=[row[1]['Chromosome'],row[1]['Strand'],int(row[1]['Start']) ,int(row[1]['End'])]

#Read Coordinates excel file, same as above but this time keep only the flam-like rows
coordinates_flamlike=pd.read_excel('data/TableS1.xlsx')
coordinates_flamlike=coordinates_flamlike[coordinates_flamlike['Region'].isin(['flamlike1','flamlike2','flamlike3','flamlike4','flamlike5','flamlike5-5p'])].reset_index(drop=True)
coordinates_flamlike['Strand']=np.where(coordinates_flamlike['Strand']=='+', 'plus', 'minus')

#Create a Dictionary, where key is the species with the build id and the values are the coordinates for the flam-like
build2coords_flamlike={}
species2build_flamlike={}
for row in coordinates_flamlike.iterrows():
    if row[1]['Species'] not in species2build_flamlike.keys():
        species2build_flamlike[row[1]['Species']]=[row[1]['Assembly']]
    else:
        update_build_list=species2build_flamlike[row[1]['Species']]+[row[1]['Assembly']]
        species2build_flamlike[row[1]['Species']]=update_build_list


for row in coordinates_flamlike.iterrows():
    build2coords_flamlike[row[1]['Species']+'_build_'+row[1]['Assembly']]=[row[1]['Chromosome'],row[1]['Strand'],int(row[1]['Start']) ,int(row[1]['End'])]

#Read Centromere regions
centromere_coordinates=pd.read_csv('data/centromere_coordinates.txt', sep='\t', index_col=0)
centromere_coordinates


#Functions
def read_proces_05files(filename,edta_directory):
    '''
    input: 1)filename - do not include "plus.bed" or "minus.bed" in the name(directory_filename)
           2)directory of the files
    output: 1)file with merged plus and minus strand information and with updated column names
    '''

    for file in os.listdir(edta_directory+'plus/'):
        if file.startswith(filename):

            plus=pd.read_csv(edta_directory+'plus/'+filename+'.bed', sep='\t', header=None)
            minus=pd.read_csv(edta_directory+'minus/'+filename+'.bed', sep='\t', header=None)

            plus.columns=['chr', 'bin_start', 'bin_end', 'transposon_overlap', 'base_overlap','bin_length', 'coverage']
            minus.columns=['chr', 'bin_start', 'bin_end', 'transposon_overlap', 'base_overlap','bin_length', 'coverage']

            plus['strand']='plus'
            minus['strand']='minus'

            plus['coverage']=plus['coverage'].astype(float)
            minus['coverage']=minus['coverage'].astype(float)

            merged=pd.concat([plus, minus])
            merged=pd.merge(left=plus,
                            right=minus,
                            left_on=['chr', 'bin_start', 'bin_end','bin_length'],
                            right_on=['chr', 'bin_start', 'bin_end','bin_length'],
                            suffixes=['_plus','_minus'])
            #X=np.array([[i] for i  in pd.concat([plus['y'],minus['y']]).tolist()])

            return(merged)


        
        
def genome_cluster_annotation(species_build_id, df):
    '''
    input:
    - Species and build id - this should match any key in the build2coords dictionary
    - The training data dataframe

    output:
    - Training data dataframe with the annotation of flam clusters
    '''
    cluster=[]
    #'Dmel_build_dm6'
    query_range=range(build2coords[species_build_id][2],build2coords[species_build_id][3])
#    for row in tqdm(df.iterrows(), total=df.shape[0]):
    for row in df.iterrows():
        if row[1]['chr']!= build2coords[species_build_id][0]:
            cluster.append(0)
        else:
#            if row[1]['bin_start'] and row[1]['bin_end'] not in query_range:
            if row[1]['bin_start'] not in query_range and row[1]['bin_end'] not in query_range:
                cluster.append(0)

            else:
                binStart=row[1]['bin_start']
                binEnd=row[1]['bin_end']
                binRange=range(binStart, binEnd)
                sbinRange=set(binRange)
                intersection_bin_cluster=sbinRange.intersection(query_range)
                if len(intersection_bin_cluster)==0:
                    print("This should not happen?")
                    cluster.append(0)
                else:
                    cluster.append(1)

    df['cluster']=cluster

    return (df)

def genome_cluster_annotation_flamlike(species_build_id, df):
    '''
    input:
    - Species and build id - this should match any key in the build2coords dictionary
    - The training data dataframe

    output:
    - Training data dataframe with the annotation of flam/flamlike clusters
    '''
    cluster=[]
    #'Dmel_build_dm6'
    query_range=range(build2coords_flamlike[species_build_id][2],build2coords_flamlike[species_build_id][3])
#    for row in tqdm(df.iterrows(), total=df.shape[0]):
    for row in df.iterrows():
        if row[1]['chr']!= build2coords_flamlike[species_build_id][0]:
            cluster.append(0)
        else:
#            if row[1]['bin_start'] and row[1]['bin_end'] not in query_range:
            if row[1]['bin_start'] not in query_range and row[1]['bin_end'] not in query_range:
                cluster.append(0)

            else:
                binStart=row[1]['bin_start']
                binEnd=row[1]['bin_end']
                binRange=range(binStart, binEnd)
                sbinRange=set(binRange)
                intersection_bin_cluster=sbinRange.intersection(query_range)
                if len(intersection_bin_cluster)==0:
                    print("This should not happen?")
                    cluster.append(0)
                else:
                    cluster.append(1)

    df['cluster_flamlike']=cluster

    return (df)

def calculate_transmat(df, number, strand, pseudo=False):
#    starprob=np.array([0.34, 0.33, 0.33])

    NC_NC=0
    NC_C=0
    NC_Cent=0
    C_NC=0
    C_C=0
    C_Cent=0
    Cent_NC=0
    Cent_C=0
    Cent_Cent=0

    NC_start = 0
    C_start = 0
    Cent_start = 0

    # Split per chromosome or not
    for chromosome in df['chr'].unique():
        df_chr=df[(df['chr']==chromosome)]

        # Calculate starting probabilities
        ends = [np.array(df_chr['region_binary'])[0],np.array(df_chr['region_binary'])[df_chr.shape[0]-1]]
        for e in ends:
            if e==0:
                NC_start+=1
            elif e==1:
                C_start+=1
            elif e==2:
                Cent_start+=1
            else:
                print("Error")

        # Calculate transition probabilities
        if strand=='plus':
            transitions_strand_one=np.vstack([df_chr['region_binary'][:-1],df_chr['region_binary'][1:]]).T
            transitions_strand_two=np.vstack([df_chr['region_binary'][:-1].replace(1,0),df_chr['region_binary'][1:].replace(1,0)]).T
        elif strand=='minus':
            transitions_strand_one=np.vstack([df_chr['region_binary'][:-1].replace(1,0),df_chr['region_binary'][1:].replace(1,0)]).T
            transitions_strand_two=np.vstack([df_chr['region_binary'][:-1],df_chr['region_binary'][1:]]).T
        else:
            print("Error: Strand has to be 'plus' or 'minus'")

        transitions=np.vstack([transitions_strand_one, np.flip(transitions_strand_two)])

        for transition in transitions:
            if transition.tolist()==[0,0]:
                NC_NC+=1
            elif transition.tolist()==[0,1]:
                NC_C+=1
            elif transition.tolist()==[0,2]:
                NC_Cent+=1
            elif transition.tolist()==[1,0]:
                C_NC+=1
            elif transition.tolist()==[1,1]:
                C_C+=1
            elif transition.tolist()==[1,2]:
                C_Cent+=1
            elif transition.tolist()==[2,0]:
                Cent_NC+=1
            elif transition.tolist()==[2,1]:
                Cent_C+=1
            elif transition.tolist()==[2,2]:
                Cent_Cent+=1
            else:
                print('error')

#    print([NC_start,C_start,Cent_start])
        
    start = NC_start+C_start+Cent_start

    if pseudo==False:
        starprob=np.array([NC_start/start, C_start/start, Cent_start/start])
        transmat=np.array([[NC_NC/(NC_NC +NC_C +NC_Cent), NC_C/(NC_NC +NC_C +NC_Cent), NC_Cent/(NC_NC +NC_C +NC_Cent)],
              [C_NC/(C_NC +C_C +C_Cent), C_C/(C_NC +C_C +C_Cent), C_Cent/(C_NC +C_C +C_Cent)],
              [Cent_NC/(Cent_NC +Cent_C +Cent_Cent),Cent_C/(Cent_NC +Cent_C +Cent_Cent), Cent_Cent/(Cent_NC +Cent_C +Cent_Cent)]])
        
    elif pseudo==True:
        NC_NC+=1e-3
        NC_C+=1e-3
        NC_Cent+=1e-3
        C_NC+=1e-3
        C_C+=1e-3
        C_Cent+=1e-3
        Cent_NC+=1e-3
        Cent_C+=1e-3
        Cent_Cent+=1e-3

        start += 3
        starprob=np.array([(NC_start+1)/start, (C_start+1)/start, (Cent_start+1)/start])
            
        transmat=np.array([[NC_NC/(NC_NC +NC_C +NC_Cent), NC_C/(NC_NC +NC_C +NC_Cent), NC_Cent/(NC_NC +NC_C +NC_Cent)],
              [C_NC/(C_NC +C_C +C_Cent), C_C/(C_NC +C_C +C_Cent), C_Cent/(C_NC +C_C +C_Cent)],
              [Cent_NC/(Cent_NC +Cent_C +Cent_Cent),Cent_C/(Cent_NC +Cent_C +Cent_Cent), Cent_Cent/(Cent_NC +Cent_C +Cent_Cent)]])
    else:
        print('Error pseudo should be either =False or =True')
    
    return(starprob, transmat)


def calculate_3emissions(df, threshold):
    emission_new=[]
    for row in df.iterrows():
        if row[1]['coverage_minus']-row[1]['coverage_plus']>threshold: # Minus strand has TE bias above threshold
            emission_new.append(1)
        elif row[1]['coverage_plus']+row[1]['coverage_minus']>threshold and row[1]['coverage_plus']-row[1]['coverage_minus']<threshold: # Both strands together are above threshold
            emission_new.append(2)
        else:
            emission_new.append(0)
    for row in df.iterrows():
        if row[1]['coverage_plus']-row[1]['coverage_minus']>threshold: # Minus strand has TE bias above threshold
            emission_new.append(1)
        elif row[1]['coverage_plus']+row[1]['coverage_minus']>threshold and row[1]['coverage_minus']-row[1]['coverage_plus']<threshold:
            emission_new.append(2)
        else:
            emission_new.append(0)
            
    return(emission_new)


def calculate_emissionprob_em3_st3(df, emission_new, strand, pseudo=False):
    if strand=='plus':
        emission_prob=pd.DataFrame(pd.concat([df['region_binary'], pd.Series([i for i in df['region_binary']]).replace(1,0)], axis=0))        
    elif strand=='minus':
        emission_prob=pd.DataFrame(pd.concat([pd.Series([i for i in df['region_binary']]).replace(1,0), df['region_binary']], axis=0))
    elif strand=='both':
        plus_strand=[]
        minus_strand=[]
        for species in df.Data.unique():
            if df[df['Data']==species]['strand'].unique()[0]=='plus':
                minus_strand=minus_strand+pd.Series([i for i in df[df['Data']==species]['region_binary']]).replace(1,0).tolist()
                plus_strand=plus_strand+df[df['Data']==species]['region_binary'].tolist()
                
            elif df[df['Data']==species]['strand'].unique()[0]=='minus':
                minus_strand=minus_strand+df[df['Data']==species]['region_binary'].tolist()
                plus_strand=plus_strand+pd.Series([i for i in df[df['Data']==species]['region_binary']]).replace(1,0).tolist()
            
        emission_prob=pd.DataFrame(pd.concat([pd.DataFrame(plus_strand), pd.DataFrame(minus_strand)], axis=0))
        
    emission_prob=emission_prob.reset_index(drop=True)
    emission_prob=pd.concat([emission_prob, pd.DataFrame(emission_new)], axis=1)
    emission_prob.columns=['region_binary','emission_new']
    emission_prob['value']=1
    emission_prob=emission_prob.groupby(['region_binary','emission_new']).count().reset_index()
    emission_prob
    try: No_cluster_0=emission_prob[(emission_prob['region_binary']==0)&(emission_prob['emission_new']==0)]['value'].values[0]
    except: No_cluster_0=0
    try: No_cluster_1=emission_prob[(emission_prob['region_binary']==0)&(emission_prob['emission_new']==1)]['value'].values[0]
    except: No_cluster_1=0
    try: No_cluster_2=emission_prob[(emission_prob['region_binary']==0)&(emission_prob['emission_new']==2)]['value'].values[0]
    except: No_cluster_2=0
    try: cluster_0=emission_prob[(emission_prob['region_binary']==1)&(emission_prob['emission_new']==0)]['value'].values[0]
    except: cluster_0=0
    try: cluster_1=emission_prob[(emission_prob['region_binary']==1)&(emission_prob['emission_new']==1)]['value'].values[0]
    except: cluster_1=0
    try: cluster_2=emission_prob[(emission_prob['region_binary']==1)&(emission_prob['emission_new']==2)]['value'].values[0]
    except: cluster_2=0
    try: centro_0=emission_prob[(emission_prob['region_binary']==2)&(emission_prob['emission_new']==0)]['value'].values[0]
    except: centro_0=0
    try: centro_1=emission_prob[(emission_prob['region_binary']==2)&(emission_prob['emission_new']==1)]['value'].values[0]
    except: centro_1=0
    try: centro_2=emission_prob[(emission_prob['region_binary']==2)&(emission_prob['emission_new']==2)]['value'].values[0]
    except: centro_2=0
    
    if pseudo==False:
        emissionmat=np.array([[No_cluster_0/(No_cluster_0+No_cluster_1+No_cluster_2), No_cluster_1/(No_cluster_0+No_cluster_1+No_cluster_2),No_cluster_2/(No_cluster_0+No_cluster_1+No_cluster_2)],
                              [cluster_0/(cluster_0+cluster_1+cluster_2), cluster_1/(cluster_0+cluster_1+cluster_2), cluster_2/(cluster_0+cluster_1+cluster_2)],
                              [centro_0/(centro_0+centro_1+centro_2), centro_1/(centro_0+centro_1+centro_2), centro_2/(centro_0+centro_1+centro_2)]])
    elif pseudo==True:
        No_cluster_0+=1
        No_cluster_1+=1
        No_cluster_2+=1
        cluster_0+=1
        cluster_1+=1
        cluster_2+=1
        centro_0+=1
        centro_1+=1
        centro_2+=1
        emissionmat=np.array([[No_cluster_0/(No_cluster_0+No_cluster_1+No_cluster_2), No_cluster_1/(No_cluster_0+No_cluster_1+No_cluster_2),No_cluster_2/(No_cluster_0+No_cluster_1+No_cluster_2)],
                              [cluster_0/(cluster_0+cluster_1+cluster_2), cluster_1/(cluster_0+cluster_1+cluster_2), cluster_2/(cluster_0+cluster_1+cluster_2)],
                              [centro_0/(centro_0+centro_1+centro_2), centro_1/(centro_0+centro_1+centro_2), centro_2/(centro_0+centro_1+centro_2)]])
                          
    return(emissionmat)
    

def train_model(number_emissions, number_states, starprob, transmat, emissionmat, all_data, emission_new, all_data_test, emission_new_test, threshold):
    model=hmm.MultinomialHMM(n_components=number_states,n_iter=100000,random_state=13)
    model.n_features=number_emissions
    model.startprob_=starprob
    model.transmat_=transmat
    model.emissionprob_=emissionmat
    
    X_train = np.atleast_2d(emission_new).T
    
    
    minus_df=all_data_test[['chr','bin_start', 'bin_end']]
    minus_df.insert(minus_df.shape[1],'strand', 'minus')
    plus_df=all_data_test[['chr','bin_start', 'bin_end']]
    plus_df.insert(plus_df.shape[1],'strand', 'plus')
    new_index=pd.concat([plus_df, minus_df]).reset_index(drop=True)
    X_test=pd.concat([new_index,pd.DataFrame(emission_new_test)], axis=1)
    
    predictions=[]
    probabilities_NC=[]
    probabilities_C=[]
    probabilities_Cent=[]
    
    for strand in ['plus','minus']:
        for chromosome in X_test['chr'].unique():
            X_test_chr=X_test[(X_test['chr']==chromosome) & (X_test['strand']==strand)]
            X_test_chr=np.atleast_2d(X_test_chr[0].tolist()).T
            predictions=predictions+ model.predict(X_test_chr).tolist()
            probabilities=model.predict_proba(X_test_chr)
            probabilities_NC=probabilities_NC+probabilities.T[0].tolist()
            probabilities_C=probabilities_C+probabilities.T[1].tolist()
            probabilities_Cent=probabilities_Cent+probabilities.T[2].tolist()
        
    X_test['pred_'+str(threshold)]=predictions
    X_test['proba_NoCluster_'+str(threshold)]=probabilities_NC
    X_test['proba_Cluster_'+str(threshold)]=probabilities_C
    X_test['proba_Centromere_'+str(threshold)]=probabilities_Cent
    X_test=X_test.rename(columns={0:'emission_'+str(threshold)})
        
    
    results_concat=X_test.set_index(['chr','bin_start', 'bin_end','strand'])
    return(results_concat)
    
