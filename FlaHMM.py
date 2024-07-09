#!/usr/bin/env python

"""
FlaHMM.py - A tool for flam-like cluster predictions based on LTR/Gypsy transposon content.
"""
_version_ = "0.9.1-beta"
_licence_ = "MIT License"

import os
import sys

import pandas as pd
import numpy as np
import sklearn
import pickle

from hmmlearn import hmm
try:  # version > 0.2.7
    from hmmlearn.hmm import CategoricalHMM as MultinomialHMM
except BaseException:  # version <= 0.2.7
    from hmmlearn.hmm import MultinomialHMM

import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns
from optparse import OptionParser

sys.path.append('scripts')
from flaHMM_functions import *


def print_intro_message():
    """ Print a welcome message """

    welcome_message = f"""=========================================================================================
FlaHMM v{_version_}

A tool for flam-like cluster predictions based on LTR/Gypsy transposon content.

If you encounter any issues or have any suggestions, please visit our GitHub issues page:
https://github.com/Hannon-lab/FlaHMM/issues
========================================================================================="""

    print(welcome_message)


def main():
    """ Main function """

    # Print welcome message
    print_intro_message()

    # Retrive user-submitted options
    parser = OptionParser()
    parser.add_option('--species', dest='species', default=None, type=str, help='name of species to predict flam-like clusters')
    parser.add_option('--bins', dest='bins', default=None, type=str, help='specify the genome bins size in kb [2.5, 5, 10]')
    parser.add_option('--threshold', dest='threshold', default=None, type=float, help='LTR/Gypsy threshold [0.025, 0.05, 0.075, 0.1, 0.2]')
    try:
        (options, args) = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    # Check that all mandatory arguments are provided
    if options.bins == None or options.threshold == None or options.species == None:
        parser.print_help()
        sys.exit(0)

    # Check that pre-computed TE annotation files are available
    bed_plus = 'bins/bins_'+str(options.bins)+'k/plus/'+str(options.species)+'.bed'
    bed_minus = 'bins/bins_'+str(options.bins)+'k/minus/'+str(options.species)+'.bed'
    if not os.path.exists(bed_plus):
        raise ValueError(bed_plus+' does not exist. Please ensure that pre-processed annotations are available.')
    if not os.path.exists(bed_minus):
        raise ValueError(bed_minus+' does not exist. Please ensure that pre-processed annotations are available.')

    # Check that pre-trained model exists
    model_name = 'models_pkl/Model_bin_' + str(options.bins) + 'k_threshold_' + str(options.threshold) + '.pkl'
    if not os.path.exists(model_name):
        raise ValueError(model_name+' does not exist. Please ensure that a pre-trained model is available.')

    # Create results folder if it doesn't exist
    if not os.path.isdir('results'):
        os.mkdir('results')

    # Format test data
    print('Formatting test data....')
    all_data = read_data(options.bins, options.species)

    # Make predictions on test data and write to file
    print('Making flaHMM predictions for ' + str(options.species) + '...')
    predictions = make_predicitions(all_data, model_name, options.threshold, [options.species])
    save_file = 'results/results_' + options.species + '_Bin_' + str(options.bins) + 'k_threshold_' + str(options.threshold) + '.txt'
    predictions.reset_index().to_csv(save_file, sep='\t')
    print('\nDone. Predictions saved in results folder as '+save_file+'.')


if __name__ == "__main__":
    main()

