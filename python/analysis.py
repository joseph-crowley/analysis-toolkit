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
import json

import helpers as hf

def main():
    # Load the C++ code that will be used to process the events
    root.gSystem.Load('../cpp/loopers/analyze_bjets.so')
    
    # Load the samples
    samples_json = '../config/samples_MC_Run2.json'
    with open(samples_json,'r') as f:
        samples = json.load(f)
    
    # Create a dictionary of TChains for each sample
    map_sample_to_category = hf.map_sample_to_category
    ch = {}

    # Loop over the samples and add the files to the TChains
    for name,sample in samples:
        for period in sample: 
            # Create a string that will be used to identify the sample
            sample_str = period + "_" + map_sample_to_category[name]

            # Create a TChain for the sample if it doesn't already exist
            if sample_str not in ch:
                ch[sample_str] = root.TChain("Events")

            # Add the files to the TChain
            for file_ in hf.get_files(sample[period]['paths']): 
                ch[sample_str].Add(file_)
    
    # Process the TChains
    for sample_str in ch:
        root.process_chain(ch[sample_str], sample_str)

if __name__ == "__main__":
    main()