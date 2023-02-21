import os
import uproot
import pandas as pd
import numpy as np
import json
import glob

samples_csv = '/home/users/crowley/tttt/CMSSW_10_6_26/src/tttt/test/samples_MC_Run2.csv'
samples_json = '../config/samples_MC_Run2.json'

data_dir = '/ceph/cms/store/group/tttt/Worker/crowley/output/Analysis_TTJetRadiation/230214_tt_bkg_Cutbased'

#samples_csv = 'samples.csv'
#samples_json = samples_csv.replace('csv','json')

# Functions to calculate the scale1fb and other metadata for a set of files and cross section
def calc_scale1fb(xs, sum_weights):
    if xs <= 0:
        return -1
    else:
        return (xs * 1000.) / sum_weights

def calculate_metadata(files, xs, debug):
    # for tttt we don't use this
    return {"xs": xs}

    # other analysis later
    nEvents = 0
    sumWeights = 0

    files_to_use = get_files(files)

    for file_ in files_to_use:
        root_file = uproot.open(file_)
        runs = root_file["Events"]
        if not runs.keys():
            results = {
                "xs" : xs,
                "scale1fb" : 1,
                "n_events" : 0,
                "sumWeights" : 0
            }
            return results
        nEntries = runs.num_entries

        nEvents += nEntries 
        # sumWeights += TODO

    # scale1fb is already in the weights, not used for this analysis
    scale1fb = calc_scale1fb(xs, sumWeights)

    results = {
        "xs" : xs,
        "scale1fb" : scale1fb,
        "n_events" : nEvents,
    }
    return results

# Function to get the paths for a given dataset and period
def get_paths(short_name, period):
    if period == '2016':
        prefix = '2016_APV' if 'APV' in short_name else '2016_NonAPV'
    else:
        prefix = period

    pattern = f'{short_name}_*_of_*.root'
    path = os.path.join(data_dir, prefix, pattern)

    return [path]

def get_files(files):
    file_list = []
    for file_ in files:
        try:
            file_list.extend(glob.glob(file_))
        except (FileNotFoundError, TypeError, OSError) as e:
            # Handle the error here, for example by printing an error message or raising a custom exception.
            print(f'Error: {e}')
    return file_list

# Load CSV file into a pandas DataFrame
df = pd.read_csv(samples_csv)

# Initialize an empty dictionary to hold the output
samples_dict = {}

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    # Extract the relevant information from the row
    dataset = row['dataset']
    short_name = row['short_name']

    if "#" in dataset:
        continue

    if short_name not in samples_dict:
        samples_dict[short_name] = {}

    period = row['period']
    cross_section = row['xsec']
    branching_ratio = row['BR']
    extra_arguments = row['extra_arguments']

    # Create a new dictionary for the period if it doesn't exist yet
    if period not in samples_dict[short_name]:
        samples_dict[short_name][period] = {'paths':[], 'metadata':{}}

    # Add the paths and metadata to the period dictionary
    file_paths = get_paths(short_name, period)
    metadata = calculate_metadata(file_paths, cross_section, 0)
    samples_dict[short_name][period]['paths'].append(file_paths)
    samples_dict[short_name][period]['metadata'].update(metadata)

# Convert the dictionary to a JSON string and save it to a file
with open(samples_json, 'w') as f:
    json.dump(samples_dict, f, indent=4)
