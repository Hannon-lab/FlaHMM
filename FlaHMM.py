import os
import sys

import pandas as pd
import numpy as np
import sklearn
from sklearn.model_selection import KFold
import pickle
from hmmlearn import hmm
try: # version > 0.2.7
	from hmmlearn.hmm import CategoricalHMM as MultinomialHMM
except: # version <= 0.2.7
	from hmmlearn.hmm import MultinomialHMM
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns
from optparse import OptionParser

sys.path.append('scripts')
from flaHMM_functions import *

def make_predicitions(all_data, model_name, threshold):
    X_test_all=pd.DataFrame()
    for species_test in species_build:
        
        X_test_all_species_test=pd.DataFrame()
        for threshold in [threshold]:
            threshold_num='threshold_'+str(threshold)
            
            #Read Test Data
            all_data=all_data_species[all_data_species['Data']==species_test]

            # Add cluster coordinates if available to help evaluation of model performance
            try: 
                strand_ref=build2coords[species_test.replace('.','_build_')][1]
            except:
                try: strand_ref=build2coords_flamlike[species_test.replace('.','_build_')][1]
                except: strand_ref=None

            #Calculate emission for the test set
            emission_new=calculate_3emissions(all_data, threshold)
            
            
            with open('models_pkl/'+model_name, "rb") as file: 
                model=pickle.load(file)
            
            
            minus_df=all_data[['chr','bin_start', 'bin_end']]
            minus_df.insert(minus_df.shape[1],'strand', 'minus')
            plus_df=all_data[['chr','bin_start', 'bin_end']]
            plus_df.insert(plus_df.shape[1],'strand', 'plus')
            new_index=pd.concat([plus_df, minus_df]).reset_index(drop=True)
            X_test=pd.concat([new_index,pd.DataFrame(emission_new)], axis=1)    
            
            
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
                    elif (strand == "minus"):
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
                    else:
                        print("Error: Unknown strand")
                
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



parser=OptionParser()
parser.add_option('--species', dest='species', default=None, type=str, help='Name of the species to predict flam,flam-like clusters, e.g. Dmel.dm6')
parser.add_option('--bins', dest='bins', default=None, type=str, help='Specify the genome bins size')
parser.add_option('--threshold', dest='threshold', default=None, type=float, help='Threshold')

(options, args)=parser.parse_args()

def introMessage():
    print('=======================================================================================================')
    print('                                              FlaHMM                                                  ')
    print('                                                                                                      ')
    print(' For queries:  Please open an issue in the github page.')
    print('=====================================================================================================\n')
    return()


species_build=[options.species]

introMessage()

# Prepare input data for test species
for bin_size in [options.bins]:
    all_data=pd.DataFrame()
    for species_train in species_build:
        print(species_train)
        species_train_build='_build_'.join(species_train.split('.'))
        all_data=read_proces_05files(species_train, 'bins/bins_'+bin_size+'k/')
    
        try:
            all_data=genome_cluster_annotation(species_train_build,all_data)
        except:
            all_data['cluster']=0
       
        all_data = all_data[['chr', 'bin_start','bin_end','cluster','coverage_plus','coverage_minus']] 
        all_data['Data']=species_train
    
        centromere_2_bins=pd.DataFrame()
        for chromosome in ['chr2L','chr2R']:
            get_centromere_coordinates=centromere_coordinates[(centromere_coordinates['species']==species_train)&(centromere_coordinates['chr']==chromosome)]
            for row in get_centromere_coordinates.iterrows():
                centromere_2_bins=pd.concat([centromere_2_bins, all_data[(all_data['chr']==chromosome)&(all_data['bin_start']>=row[1]['bin_start'])&(all_data['bin_end']<=row[1]['bin_end'])]])
            
        centromere_3_bins=pd.DataFrame()
        for chromosome in ['chr3L','chr3R']:
            get_centromere_coordinates=centromere_coordinates[(centromere_coordinates['species']==species_train)&(centromere_coordinates['chr']==chromosome)]
            for row in get_centromere_coordinates.iterrows():
                centromere_3_bins=pd.concat([centromere_3_bins, all_data[(all_data['chr']==chromosome)&(all_data['bin_start']>=row[1]['bin_start'])&(all_data['bin_end']<=row[1]['bin_end'])]])
          
        flam_bins=all_data[all_data['cluster']==1]
#        print(flam_bins.shape)
    
        region=[]
        region_binary=[]
        for i in tqdm(all_data.index.tolist()):
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
            
        all_data['region']=region
        all_data['region_binary']=region_binary
        all_data_species = all_data
       
 
# Create results folder if it doesn't exist    
if not os.path.isdir('results'):
    os.mkdir('results')
    
for threshold in [options.threshold]:
    print('Making flaHMM predictions for '+str(species_build))
    get_bin_size=options.bins+'k'
    model_make_pred='Model_bin_'+get_bin_size+'_threshold_'+str(threshold)+'.pkl'
    results_predictions=make_predicitions(all_data,model_make_pred, threshold)
    results_predictions.reset_index().to_csv('results/results_'+options.species+'_Bin_'+get_bin_size+'_threshold_'+str(threshold)+'.txt', sep='\t')
    print('Predictions saved in results folder!')
    
