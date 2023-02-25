import json
import glob

def load_sample_map(tag):
    with open(f'../config/sample_map_{tag}.json','r') as f:
        sample_map = json.load(f)
    sample_map_inv = {}
    for k,l in sample_map.items():
        for s in l:
            sample_map_inv[s] = k
    return sample_map, sample_map_inv

#map_category_to_sample, map_sample_to_category = load_sample_map()

def get_samples_by_category(categories, files):
    samples = {}
    for category in categories:
        samples.update({category:[]})
        for f in files:
            if "Data" in category:
                base_sample = f.split('/')[-1].split('_')[0]
            else:
                base_sample = "_".join(f.split('/')[-1].replace("_ext","").split("_of")[0].split("_")[:-1])
            cat = get_category(base_sample, map_sample_to_category)
            if cat == category:
                samples[cat].append(f)

        if not samples[category]:
            print(f'no samples for {category}')
            del samples[category]
    return samples


def get_files(files_with_wildcard):
    file_list = []
    for file_ in files_with_wildcard:
        try:
            file_list.extend(glob.glob(file_))
        except (FileNotFoundError, TypeError, OSError) as err:
            print(f'Error retrieving files for {file_}: {err.args[0]}')
    return file_list

def get_category(base_sample, map_sample_to_category):
    if 'Run20' in base_sample: 
        return "Data_" + base_sample[3:-1]
    if base_sample in map_sample_to_category.keys():
        return map_sample_to_category[base_sample]
    return 'Ignore'
