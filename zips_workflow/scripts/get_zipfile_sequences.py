from zipfile import ZipFile
from glob import glob
import os
import json


all_sequences = set()

for f in glob(os.path.join(snakemake.input.zip_dir,'*')):

    with ZipFile(f,'r') as zipfile:

        file_list = zipfile.namelist()
    
        print(f'parsing {len(file_list)} files in the zip')

        #get set of subjects per zipfile
        for fname in file_list:
            split_fname = fname.split('/')

            sequence= split_fname[2]

            all_sequences.add(sequence)

#categorize with simple heuristics:
mapping=dict()

for seq in all_sequences:
    seql=seq.lower() #for case-insensitive
    if ('t1' in seql) or ('mprage' in seql) or ('spgr' in seql) or ('bravo' in seql) or ('mpr' in seql and 'sag' in seql):
        mapping[seq]='T1w'
    elif ('dti' in seql) or ('diff' in seql) or ('gated' in seql) or ('dwi' in seql):
        mapping[seq]='dwi'
    elif 'flair' in seql:
        mapping[seq]='FLAIR'
    elif ('rs' in seql) or ('resting' in seql):
        mapping[seq]='bold'
    elif ('nm' in seql) or ('mt' in seql):
        mapping[seq]='MT'
    elif ('pd' in seql) or ('dual' in seql):
        mapping[seq]='PDT2'
    elif 't2' in seql:
        mapping[seq]='T2w'
    elif 'swan' in seql:
        mapping[seq]='SWI'
    else:
        mapping[seq]='UNKNOWN'
        

#save to json
all_sequences = sorted(list(all_sequences))


out_dict = dict()
out_dict['all_sequences'] = all_sequences
out_dict['auto_map'] = mapping

#save both mappings to the json file
with open(snakemake.output[0], "w") as outfile:
    json.dump(out_dict, outfile, indent=4)

