from zipfile import ZipFile
from glob import glob
import os
import json

zipfiles_to_subjects = dict()

all_subjects = set()

for f in glob(os.path.join(snakemake.input.zip_dir,'*')):

    with ZipFile(f,'r') as zipfile:

        file_list = zipfile.namelist()
    
        print(f'parsing {len(file_list)} files in the zip')

        #get set of subjects per zipfile
        subjects = set()
        for fname in file_list:
            split_fname = fname.split('/')

            subject = split_fname[1]

            subjects.add(subject)

            all_subjects.add(subject)

        zipfiles_to_subjects[f] = sorted(list(subjects))
        



all_subjects = sorted(list(all_subjects))

subjects_to_zipfiles = dict()

#now we have the list of subjects per zipfile
#create the opposite mapping, list of zipfiles per subject
for subj in all_subjects:
    subjects_to_zipfiles[subj] = list()
    for zipfile,subjects in zipfiles_to_subjects.items():
        if subj in subjects:
            subjects_to_zipfiles[subj].append(zipfile)


out_dict = dict()
out_dict['subjects_to_zipfiles'] = subjects_to_zipfiles
out_dict['zipfiles_to_subjects'] = zipfiles_to_subjects
out_dict['all_subjects'] = all_subjects

#save both mappings to the json file
with open(snakemake.output[0], "w") as outfile:
    json.dump(out_dict, outfile, indent=4)

