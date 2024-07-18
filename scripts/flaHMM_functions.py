""" Help functions for FlaHMM training/use """

import os

import sklearn
import pandas as pd
import numpy as np
import pickle

from hmmlearn import hmm
try:  # version > 0.2.7
    from hmmlearn.hmm import CategoricalHMM as MultinomialHMM
except:  # version <= 0.2.7
    from hmmlearn.hmm import MultinomialHMM

import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns


def get_flam_dictionary():
    """
    Reads data/TableS1.xlsx and returns build2coords
    
    output:
    - build2coord - Dictionary with species_build -> flam coordinates as chromosome, strand, start, end
    """

    # Read coordinates of flam-syntenic clusters
    coordinates=pd.read_excel('data/TableS1.xlsx')
    coordinates=coordinates[coordinates['Region'].isin(['flam','flam-5p'])].reset_index(drop=True) # Prefer 5' end over 3' end
    coordinates['Strand']=np.where(coordinates['Strand']=='+', 'plus', 'minus')

    # Dictionary with species_assembly -> flam coordinate
    build2coords={}
    for row in coordinates.iterrows():
        build2coords[row[1]['Species']+'_build_'+row[1]['Assembly']]=[row[1]['Chromosome'],row[1]['Strand'],int(row[1]['Start']) ,int(row[1]['End'])]

    return(build2coords)


def get_flamlike_dictionary():
    """
    Reads data/TableS1.xlsx and returns build2coords_flamlike
    
    output:
    - build2coord_flamlike - Dictionary with species_build -> flamlike coordinates as chromosome, strand, start, end
    """

    # Read coordinates of flam-syntenic clusters
    coordinates_flamlike=pd.read_excel('data/TableS1.xlsx')
    coordinates_flamlike=coordinates_flamlike[coordinates_flamlike['Region'].isin(['flamlike1','flamlike2','flamlike3','flamlike4','flamlike5','flamlike5-5p'])].reset_index(drop=True)
    coordinates_flamlike['Strand']=np.where(coordinates_flamlike['Strand']=='+', 'plus', 'minus')

    # Dictionary with species_assembly -> flamlike coordinate
    build2coords_flamlike={}
    for row in coordinates_flamlike.iterrows():
        build2coords_flamlike[row[1]['Species']+'_build_'+row[1]['Assembly']]=[row[1]['Chromosome'],row[1]['Strand'],int(row[1]['Start']) ,int(row[1]['End'])]

    return(build2coords_flamlike)


def read_centromere_coordinates():
    """
    Read data/centromere_coordinates.txt and return centromere_coordinates

    output:
    - centromere_coordinates - Pandas object with centromere coordinates (species, chr, start, end, bin_start, bin_end, length)
    """

    centromere_coordinates=pd.read_csv('data/centromere_coordinates.txt', sep='\t', index_col=0)
    centromere_coordinates
    
    return(centromere_coordinates)


def read_proces_05files(filename, directory):
    """
    Read BED-like files with TE coverage per genomic bin

    input: 
    - filename - sample name excpluding the .bed" part
    - directory - directory of the files (excluding the plus/ and minus/ part)

    output: 
    - merged - Pandas object with columns: chr, bin_start, bin_end, bin_length, 
               transposon_overlap_plus, base_overlap_plus, coverage_plus, strand_plus, 
               transposon_overlap_minus, base_overlap_minus, coverage_minus, strand_minus
    """

    # Read TEs on the plus and minus strand separately
    plus=pd.read_csv(directory+'plus/'+filename+'.bed', sep='\t', header=None)
    minus=pd.read_csv(directory+'minus/'+filename+'.bed', sep='\t', header=None)

    plus.columns=['chr', 'bin_start', 'bin_end', 'transposon_overlap', 'base_overlap','bin_length', 'coverage']
    minus.columns=['chr', 'bin_start', 'bin_end', 'transposon_overlap', 'base_overlap','bin_length', 'coverage']

    plus['strand']='plus'
    minus['strand']='minus'

    plus['coverage']=plus['coverage'].astype(float)
    minus['coverage']=minus['coverage'].astype(float)

    merged=pd.merge(left=plus,
                    right=minus,
                    left_on=['chr', 'bin_start', 'bin_end','bin_length'],
                    right_on=['chr', 'bin_start', 'bin_end','bin_length'],
                    suffixes=['_plus','_minus'])

    return(merged)
        
        
def genome_cluster_annotation(species_build_id, df):
    """
    Takes a pandas DataFrame and adds flam annotations as 'cluster'

    input:
    - species_build_id - Name of assembly; this should match any key in the build2coords dictionary
    - df - Pandas DataFrame

    output:
    - df - Pandas DataFrame with 'cluster' column added
    """

    build2coords = get_flam_dictionary()

    cluster=[]
    query_range=range(build2coords[species_build_id][2],build2coords[species_build_id][3])
    for row in df.iterrows():
        if row[1]['chr'] != build2coords[species_build_id][0]:
            cluster.append(0)
        else:
            if row[1]['bin_start'] not in query_range and row[1]['bin_end'] not in query_range:
                cluster.append(0)
            else:
                cluster.append(1)

    df['cluster']=cluster
    return (df)


def genome_cluster_annotation_flamlike(species_build_id, df):
    """
    Takes a pandas DataFrame and adds flam annotations as 'cluster'

    input:
    - species_build_id - Name of assembly; this should match any key in the build2coords_flamlike dictionary
    - df - Pandas DataFrame

    output:
    - df - Pandas DataFrame with 'cluster_flamlike' column added
    """

    build2coords_flamlike = get_flamlike_dictionary()

    cluster=[]
    query_range=range(build2coords_flamlike[species_build_id][2],build2coords_flamlike[species_build_id][3])
    for row in df.iterrows():
        if row[1]['chr'] != build2coords_flamlike[species_build_id][0]:
            cluster.append(0)
        else:
            if row[1]['bin_start'] not in query_range and row[1]['bin_end'] not in query_range:
                cluster.append(0)
            else:
                cluster.append(1)

    df['cluster_flamlike']=cluster
    return (df)


def read_data(bins, species_build):
    """
    Read transposon coverage for selected sample and bin size.

    input:
    - bins - Bin size in kb
    - species_build - Name of genome assembly

    output:
    - df - Pandas DataFrame with the annotation of centromeres and flam clusters
    """

    # Prepare input data for test species
    df = pd.DataFrame()
    df = read_proces_05files(species_build, 'bins/bins_' + bins + 'k/')

    # Retrieve flam annotations if available
    try:
        df = genome_cluster_annotation(str(species_build)+"_build", df)
    except:
        df['cluster'] = 0

    # Keep only selected columns
    df = df[['chr', 'bin_start', 'bin_end', 'cluster', 'coverage_plus', 'coverage_minus']]
    df['Data'] = species_build

    centromere_coordinates = read_centromere_coordinates()

    # Construct lists of centromere and cluster coordinates
    centromere_2_bins = pd.DataFrame()
    for chromosome in ['chr2L', 'chr2R']:
        get_centromere_coordinates = centromere_coordinates[(
            centromere_coordinates['species'] == species_build) & (centromere_coordinates['chr'] == chromosome)]
        for row in get_centromere_coordinates.iterrows():
            centromere_2_bins = pd.concat([centromere_2_bins, df[(df['chr'] == chromosome) & (
                df['bin_start'] >= row[1]['bin_start']) & (df['bin_end'] <= row[1]['bin_end'])]])

    centromere_3_bins = pd.DataFrame()
    for chromosome in ['chr3L', 'chr3R']:
        get_centromere_coordinates = centromere_coordinates[(
            centromere_coordinates['species'] == species_build) & (centromere_coordinates['chr'] == chromosome)]
        for row in get_centromere_coordinates.iterrows():
            centromere_3_bins = pd.concat([centromere_3_bins, df[(df['chr'] == chromosome) & (
                df['bin_start'] >= row[1]['bin_start']) & (df['bin_end'] <= row[1]['bin_end'])]])

    flam_bins = df[df['cluster'] == 1]

    # Create annotations based on overlap to centromere or cluster bins
    region = []
    region_binary = []
    for i in tqdm(df.index.tolist()):
        if i in flam_bins.index.tolist():
            region.append('flam')
            region_binary.append(1)
        elif i in centromere_2_bins.index.tolist():
            region.append('centromere')
            region_binary.append(2)
        elif i in centromere_3_bins.index.tolist():
            region.append('centromere')
            region_binary.append(2)
        else:
            region.append('none')
            region_binary.append(0)

    df['region'] = region
    df['region_binary'] = region_binary
    return(df)


def make_predicitions(df, model_name, threshold, species_build):
    """
    Function to make predictions using a pre-trained model across a used provided test set

    input:
    - df - Pandas DataFrame with training data (can be one or many species)
    - model_name - Saved pre-trained model
    - threshold - Threshold to be tested
    - species_build - List of species on which to make predictions

    output:
    - X_test_all - Pandas DataFrame with predictions added
    """

    build2coords = get_flam_dictionary()
    build2coords_flamlike = get_flamlike_dictionary()

    X_test_all=pd.DataFrame() # Save all results

    for species_test in species_build:

        X_test_all_species_test=pd.DataFrame() # Save all results for a given species_build

        # Select only relevant genome assembly
        all_data=df[df['Data']==species_test]

        # Add cluster coordinates if available to help evaluation of model performance
        try:
            strand_ref=build2coords[species_test.replace('.','_build_')][1]
        except:
            try: strand_ref=build2coords_flamlike[species_test.replace('.','_build_')][1]
            except: strand_ref=None

        # Calculate emission for the test set
        emissions=calculate_3emissions(all_data, threshold)

        # Load saved model
        with open(model_name, "rb") as file:
            model=pickle.load(file)

        minus_df=all_data[['chr','bin_start', 'bin_end']]
        minus_df.insert(minus_df.shape[1],'strand', 'minus')
        plus_df=all_data[['chr','bin_start', 'bin_end']]
        plus_df.insert(plus_df.shape[1],'strand', 'plus')
        df_new_index=pd.concat([plus_df, minus_df]).reset_index(drop=True)
        X_test=pd.concat([df_new_index,pd.DataFrame(emissions)], axis=1)

        predictions=[]
        probabilities_NC=[]
        probabilities_C=[]
        probabilities_Cent=[]

        for strand in ['plus','minus']:
            for chromosome in X_test['chr'].unique():
                X_test_chr=X_test[(X_test['chr']==chromosome) & (X_test['strand']==strand)]
                X_test_chr=np.atleast_2d(X_test_chr[0].tolist()).T

                if (strand == "plus"):
                    predictions=predictions + model.predict(X_test_chr).tolist()
                    probabilities=model.predict_proba(X_test_chr)
                    probabilities_NC=probabilities_NC + probabilities.T[0].tolist()
                    probabilities_C=probabilities_C + probabilities.T[1].tolist()
                    probabilities_Cent=probabilities_Cent + probabilities.T[2].tolist()
                else:
                    # Swap and re-swap to retain orientation
                    predictions_tmp = model.predict(np.flip(X_test_chr)).tolist()
                    predictions_tmp.reverse()
                    predictions=predictions + predictions_tmp

                    probabilities_tmp = model.predict_proba(np.flip(X_test_chr)).tolist()
                    probabilities_tmp.reverse()
                    probabilities_tmp = np.array(probabilities_tmp).T.tolist()
                    probabilities_NC=probabilities_NC+probabilities_tmp[0]
                    probabilities_C=probabilities_C+probabilities_tmp[1]
                    probabilities_Cent=probabilities_Cent+probabilities_tmp[2]

        X_test['pred_'+str(threshold)]=predictions
        X_test['proba_NoCluster_'+str(threshold)]=probabilities_NC
        X_test['proba_Cluster_'+str(threshold)]=probabilities_C
        X_test['proba_Centromere_'+str(threshold)]=probabilities_Cent
        X_test=X_test.rename(columns={0:'emission_'+str(threshold)})

        X_test['species_test']=species_test

        if strand_ref=='plus':
            region_binary_test=all_data['region_binary'].tolist()+all_data['region_binary'].replace(1,0).tolist()
        elif strand_ref=='minus':
            region_binary_test=all_data['region_binary'].replace(1,0).tolist()+all_data['region_binary'].tolist()
        elif  strand_ref=='both':
            region_binary_test=all_data['region_binary'].tolist()+all_data['region_binary'].tolist()
        else:
            region_binary_test=all_data['region_binary'].tolist()+all_data['region_binary'].tolist()

        X_test['region_binary']=region_binary_test

        X_test_all_species_test=pd.concat([X_test_all_species_test,X_test.set_index(['chr', 'bin_start', 'bin_end', 'strand','species_test','region_binary'])], axis=1)

        X_test_all=pd.concat([X_test_all, X_test_all_species_test], axis=0)

    return(X_test_all)


def calculate_transmat(df, number, strand, pseudo=False):
    """
    Calculate transition and starting matrices

    input:
    - df - Pandas DataFrame
    - number - Number of states; currently assumed to be 3
    - strand - Strand of flam cluster to be used for training
    - pseudo - True/False indicating whether to add psedo counts

    output:
    - starting_probs
    - transition_matrix
    """

    # Transition matrix (from_to)
    NC_NC=0
    NC_C=0
    NC_Cent=0
    C_NC=0
    C_C=0
    C_Cent=0
    Cent_NC=0
    Cent_C=0
    Cent_Cent=0

    # Starting probabilities
    NC_start = 0
    C_start = 0
    Cent_start = 0

    # Split per chromosome
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
        

    # Add pseudocounts if selected
    if pseudo==True:
        NC_NC += 1e-3
        NC_C += 1e-3
        NC_Cent += 1e-3
        C_NC += 1e-3
        C_C += 1e-3
        C_Cent += 1e-3
        Cent_NC += 1e-3
        Cent_C += 1e-3
        Cent_Cent += 1e-3

        NC_start += 1
        C_start += 1
        Cent_start += 1
      
    # Set starting probabilities 
    start = NC_start+C_start+Cent_start
    starting_probs=np.array([NC_start/start, C_start/start, Cent_start/start])

    # Set transition matrix     
    transition_matrix=np.array([[NC_NC/(NC_NC +NC_C +NC_Cent), NC_C/(NC_NC +NC_C +NC_Cent), NC_Cent/(NC_NC +NC_C +NC_Cent)],
          [C_NC/(C_NC +C_C +C_Cent), C_C/(C_NC +C_C +C_Cent), C_Cent/(C_NC +C_C +C_Cent)],
          [Cent_NC/(Cent_NC +Cent_C +Cent_Cent),Cent_C/(Cent_NC +Cent_C +Cent_Cent), Cent_Cent/(Cent_NC +Cent_C +Cent_Cent)]])
    
    return(starting_probs, transition_matrix)


def calculate_3emissions(df, threshold):
    """
    Calculate emissions for a given sample

    input:
    - df - Pandas DataFrame
    - threshold - LTR/Gypsy threshold to be applied

    output:
    - emissions
    """

    emissions=[]
    for row in df.iterrows():
        if row[1]['coverage_minus']-row[1]['coverage_plus']>threshold: # Minus strand has TE bias above threshold
            emissions.append(1)
        elif row[1]['coverage_plus']+row[1]['coverage_minus']>threshold and row[1]['coverage_plus']-row[1]['coverage_minus']<threshold: # Both strands together are above threshold
            emissions.append(2)
        else:
            emissions.append(0)
    for row in df.iterrows():
        if row[1]['coverage_plus']-row[1]['coverage_minus']>threshold: # Minus strand has TE bias above threshold
            emissions.append(1)
        elif row[1]['coverage_plus']+row[1]['coverage_minus']>threshold and row[1]['coverage_minus']-row[1]['coverage_plus']<threshold:
            emissions.append(2)
        else:
            emissions.append(0)
            
    return(emissions)


def calculate_emissionprob_em3_st3(df, emissions, strand, pseudo=False):
    """
    Calculate emission matrix

    input:
    - df - Pandas DataFrame; must be annotated with region_binary
    - emissions - List of emissions
    - strand - Strand of flam cluster to be used for training
    - pseudo - True/False indicating whether to add psedo counts

    output:
    - emission_matrix
    """

    if strand=='plus':
        all_data=pd.DataFrame(pd.concat([df['region_binary'], pd.Series([i for i in df['region_binary']]).replace(1,0)], axis=0))        
    elif strand=='minus':
        all_data=pd.DataFrame(pd.concat([pd.Series([i for i in df['region_binary']]).replace(1,0), df['region_binary']], axis=0))
    else:
        print("Error: strand has to be plus or minus")
        
    # Add emission column to DataFrame
    all_data=all_data.reset_index(drop=True)
    all_data=pd.concat([all_data, pd.DataFrame(emissions)], axis=1)
    all_data.columns=['region_binary','emissions']

    # Add value column that will be used for counting the number of occurences of each combination 
    all_data['value']=1
    all_data=all_data.groupby(['region_binary','emissions']).count().reset_index()

    # Save counts for each state and emission combination; use except to set to 0 if combination does not exist
    try: none_0=all_data[(all_data['region_binary']==0)&(all_data['emissions']==0)]['value'].values[0]
    except: none_0=0
    try: none_1=all_data[(all_data['region_binary']==0)&(all_data['emissions']==1)]['value'].values[0]
    except: none_1=0
    try: none_2=all_data[(all_data['region_binary']==0)&(all_data['emissions']==2)]['value'].values[0]
    except: none_2=0
    try: cluster_0=all_data[(all_data['region_binary']==1)&(all_data['emissions']==0)]['value'].values[0]
    except: cluster_0=0
    try: cluster_1=all_data[(all_data['region_binary']==1)&(all_data['emissions']==1)]['value'].values[0]
    except: cluster_1=0
    try: cluster_2=all_data[(all_data['region_binary']==1)&(all_data['emissions']==2)]['value'].values[0]
    except: cluster_2=0
    try: centro_0=all_data[(all_data['region_binary']==2)&(all_data['emissions']==0)]['value'].values[0]
    except: centro_0=0
    try: centro_1=all_data[(all_data['region_binary']==2)&(all_data['emissions']==1)]['value'].values[0]
    except: centro_1=0
    try: centro_2=all_data[(all_data['region_binary']==2)&(all_data['emissions']==2)]['value'].values[0]
    except: centro_2=0
    
    if pseudo==True:
        none_0+=1
        none_1+=1
        none_2+=1
        cluster_0+=1
        cluster_1+=1
        cluster_2+=1
        centro_0+=1
        centro_1+=1
        centro_2+=1

    emission_matrix=np.array([[none_0/(none_0+none_1+none_2), none_1/(none_0+none_1+none_2),none_2/(none_0+none_1+none_2)],
                              [cluster_0/(cluster_0+cluster_1+cluster_2), cluster_1/(cluster_0+cluster_1+cluster_2), cluster_2/(cluster_0+cluster_1+cluster_2)],
                              [centro_0/(centro_0+centro_1+centro_2), centro_1/(centro_0+centro_1+centro_2), centro_2/(centro_0+centro_1+centro_2)]])
                          
    return(emission_matrix)
    

def train_model(number_emissions, number_states, starting_probs, transition_matrix, emission_matrix, df_train, emissions_train, df_test, emissions_test, threshold):
    """ 
    Defines a model with the given parameters and returns predcitions based on this model

    input:
    - number_emissions - Number of emissions (3)
    - number_states - Number of states (3)
    - starting_probs - Starting probabilities
    - transition_matrix - Transition matrix
    - emission_matrix - Emission matrix
    - df_train - Pandas DataFrame for training set (NOT USED)
    - emissions_train - Emissions for training set (NOT USED)
    - df_test - Pandas DataFrame for test set
    - emissions_test - Emissions for test set
    - threshold - Threshold to use for LTR/Gypsy content; used for column names

    output:
    - df - Pandas DataFrame for test set with predictions
    """

    # Create a CategoricalHMM model with the provided starting probabilites and transition and emission matrices
    model = hmm.MultinomialHMM(n_components=number_states, n_iter=100000, random_state=13)
    model.n_features = number_emissions
    model.startprob_= starting_probs
    model.transition_matrix_= transition_matrix
    model.emissionprob_ = emission_matrix

    # Format test data into plus and minus DataFrames
    minus_df = df_test[['chr','bin_start', 'bin_end']]
    minus_df.insert(minus_df.shape[1],'strand', 'minus')
    plus_df = df_test[['chr','bin_start', 'bin_end']]
    plus_df.insert(plus_df.shape[1],'strand', 'plus')

    # Create a new DataFrame to be returned
    df_new_index = pd.concat([plus_df, minus_df]).reset_index(drop=True)
    df = pd.concat([df_new_index, pd.DataFrame(emissions_test)], axis=1)
   
    # Lists used to save our predictions 
    predictions=[]
    probabilities_NC=[]
    probabilities_C=[]
    probabilities_Cent=[]

    # Make predictions for each strand and chromosome, one at a time, and append to lists
    for strand in ['plus','minus']:
        for chromosome in df['chr'].unique():
            df_chr=df[(df['chr']==chromosome) & (df['strand']==strand)]
            df_chr=np.atleast_2d(df_chr[0].tolist()).T
            predictions += model.predict(df_chr).tolist()
            probabilities = model.predict_proba(df_chr)
            probabilities_NC += probabilities.T[0].tolist()
            probabilities_C += probabilities.T[1].tolist()
            probabilities_Cent += probabilities.T[2].tolist()

    # Annotate final object with results 
    df['pred_'+str(threshold)]=predictions
    df['proba_NoCluster_'+str(threshold)]=probabilities_NC
    df['proba_Cluster_'+str(threshold)]=probabilities_C
    df['proba_Centromere_'+str(threshold)]=probabilities_Cent
    df = df.rename(columns={0:'emission_'+str(threshold)})
    df = df.set_index(['chr','bin_start', 'bin_end','strand'])

    return(df)
    
