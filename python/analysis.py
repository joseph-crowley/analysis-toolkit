"""
analysis.py

This script provides the main functionality for performing event-by-event analysis of particle physics data from CMS.

Functions:
    - main: This function is the main entry point for the script. It loads the samples, creates the TChains, and processes the events.
        depends on: C++ process_chain: This function takes a chain of ROOT files and processes them event-by-event.

Author: Joe Crowley
Date: 2023-02-20
"""

import ROOT as root
from ROOT import gROOT

import os
import json

import helpers as hf

TAG = 'mc'

def main():
    # Load the C++ code that will be used to process the events
    status = root.gSystem.Load('../cpp/loopers/loopers_cpp.so')

    input(f'status {status}')

    # Load the samples
    samples_json = f'../config/samples_{TAG}_Run2.json'
    with open(samples_json,'r') as f:
        samples = json.load(f)

    # Create a dictionary of TChains for each sample
    _,map_sample_to_category = hf.load_sample_map(TAG)
    categories = list(set(map_sample_to_category.values()))
    ch = {}

    # Loop over the samples and add the files to the TChains
    for name,sample in samples.items():
        if name not in map_sample_to_category:
            print(f'{name} not in sample map. Skipping...')
            continue

        category = map_sample_to_category[name]

        for period in sample:
            # Create a string that will be used to identify the sample
            sample_str = period + "_" + category

            # Create a TChain for the sample if it doesn't already exist
            if sample_str not in ch:
                ch[sample_str] = root.TChain("Events")

            # Add the files to the TChain
            files_by_sample = hf.get_files(sample[period]['paths'][0])
            
            for file_ in files_by_sample:
                ch[sample_str].Add(file_)

    # Process the TChains
    for sample_str, chain in ch.items():
        root.process_chain(chain, sample_str)

if __name__ == "__main__":
    main()
