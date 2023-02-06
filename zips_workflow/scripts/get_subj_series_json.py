from zipfile import ZipFile
from glob import glob
import os
import json
import re

all_series_descs = set()

with ZipFile(snakemake.input.zip_file,'r') as zipfile:

    file_list = zipfile.namelist()

    print(f'parsing {len(file_list)} files in the zip')

    #get set of subjects per zipfile
    for fname in file_list:
        print(fname)
        dicom_filename = fname.split('/')[-1]

        print(dicom_filename)
        #get the string for the series description
        series_desc = re.match(r'PPMI_[0-9]+_(.+)_br_.+',dicom_filename)

        if series_desc==None:
            continue
        print(series_desc)
        all_series_descs.add(series_desc.group(1))

all_series_descs = sorted(list(all_series_descs))

#categorize with simple heuristics:
mapping=dict()

for seq in all_series_descs:
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
out_dict = dict()
out_dict['all_series_descs'] = all_series_descs
out_dict['auto_map'] = mapping

#save both mappings to the json file
with open(snakemake.output[0], "w") as outfile:
    json.dump(out_dict, outfile, indent=4)

