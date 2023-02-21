"""
analysis.py

This script provides the main functions for performing event-by-event analysis of particle physics data from CMS.

Functions:
    - process_chain: This function takes a chain of ROOT files and processes them event-by-event.

Author: Joe Crowley
Date: 2023-02-20
"""

import ROOT as root
from ROOT import gROOT
import json

import helpers as hf

r.gSystem.Load('../cpp/loopers/analyze_bjets.so')

samples_json = '../config/samples_MC_Run2.json'
with open(samples_json,'r') as f:
    samples = json.load(f)

for name,sample in samples:
    ch = root.TChain("Events")

    for file_ in hf.get_files(sample['2017']['paths']): ch.Add(file_)

    root.process_chain(ch, '2017_'+name)
