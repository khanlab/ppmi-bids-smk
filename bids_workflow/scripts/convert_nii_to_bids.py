from snakemake.io import glob_wildcards
from snakemake.io import strip_wildcard_constraints
import os
import shutil
import re

#search pattern for nii.gz files in the subj_nii folder (output from dcm2niix)
input_string_with_constraints = f'{snakemake.input.nii_dir}/sub-{snakemake.wildcards.subject}_ses-{{session,[0-9]+}}_{{desc,[^.]+}}.{{ext,nii.gz}}'

#glob to get the wildcards
sessions,descs,exts = glob_wildcards(input_string_with_constraints)

#for formatting:
input_string = strip_wildcard_constraints(input_string_with_constraints)

subject = snakemake.wildcards.subject

output_string = { 
            't1': f'{snakemake.output.subj_bids_dir}/ses-{{session}}/anat/sub-{subject}_ses-{{session}}_run-{{run}}_T1w.{{ext}}',
            'dwi': f'{snakemake.output.subj_bids_dir}/ses-{{session}}/dwi/sub-{subject}_ses-{{session}}_run-{{run}}_dwi.{{ext}}',
            'bold': f'{snakemake.output.subj_bids_dir}/ses-{{session}}/func/sub-{subject}_ses-{{session}}_run-{{run}}_bold.{{ext}}',
            'flair': f'{snakemake.output.subj_bids_dir}/ses-{{session}}/anat/sub-{subject}_ses-{{session}}_run-{{run}}_FLAIR.{{ext}}',
            'nm': f'{snakemake.output.subj_bids_dir}/ses-{{session}}/anat/sub-{subject}_ses-{{session}}_run-{{run}}_MTw.{{ext}}',
            'pdt2': f'{snakemake.output.subj_bids_dir}/ses-{{session}}/anat/sub-{subject}_ses-{{session}}_echo-{{echo}}_PDT2.{{ext}}',
            'pd': f'{snakemake.output.subj_bids_dir}/ses-{{session}}/anat/sub-{subject}_ses-{{session}}_run-{{run}}_PDw.{{ext}}',
            't2': f'{snakemake.output.subj_bids_dir}/ses-{{session}}/anat/sub-{subject}_ses-{{session}}_run-{{run}}_T2w.{{ext}}',
            'swan': f'{snakemake.output.subj_bids_dir}/ses-{{session}}/anat/sub-{subject}_ses-{{session}}_run-{{run}}_T2starw.{{ext}}',
            }

sessions_set = sorted(list(set(sessions)))
print(f'sessions available for this subject: {sessions_set}')

#initialize the mapping, indexing by session as well so we can count the runs within a session
seq_dict= { seq:{ses: [] for ses in sessions_set  } for seq in output_string.keys()  }
    
for ses,desc,ext in zip(sessions,descs,exts): #extension is constrained to be nii.gz, so can leave this out

    in_file = input_string.format(session=ses,desc=desc,ext=ext)

    

    descl=desc.lower() #for case-insensitive
    if ('t1' in descl) or ('mprage' in descl) or ('spgr' in descl) or ('bravo' in descl) or ('mpr' in descl and 'sag' in descl):
        seq_dict['t1'][ses].append({'session': ses,'desc':desc})

    elif ('dti' in descl) or ('diff' in descl) or ('gated' in descl) or ('dwi' in descl):

        #check if bvec and bval exist:
        if os.path.exists(input_string.format(session=ses,desc=desc,ext='bvec')) and os.path.exists(input_string.format(session=ses,desc=desc,ext='bval')):
            seq_dict['dwi'][ses].append({'session': ses,'desc': desc})

    elif 'flair' in descl:
        seq_dict['flair'][ses].append({'session': ses,'desc':desc})
    elif ('rs' in descl) or ('resting' in descl):
        seq_dict['bold'][ses].append({'session': ses,'desc':desc})
    elif ('nm' in descl) or ('mt' in descl):
        seq_dict['nm'][ses].append({'session': ses,'desc':desc})
    elif ('pd' in descl) or ('dual' in descl):
        #check echo number if it is present:
        echo_match = re.search(r'e([0-9]+)$',desc)
        if echo_match is not None:
            echo = echo_match.group(1)
            seq_dict['pdt2'][ses].append({'session': ses,'desc':desc,'echo':echo})
        else:
            seq_dict['pd'][ses].append({'session': ses,'desc':desc})
    elif 't2' in descl:
        seq_dict['t2'][ses].append({'session': ses,'desc':desc})
    elif 'swan' in descl:
        seq_dict['swan'][ses].append({'session': ses,'desc':desc})
    else:
        print('Cannot find matching sequence for sub-{subject}, ses-{ses}: {desc}')

print(seq_dict)

#loop through types of sequences:
for seq in seq_dict.keys():
    for ses in sessions_set:
        #loop through runs:
        for run,wildcards in enumerate(seq_dict[seq][ses],start=1):

            #glob to get all the extensions for a scan:
            in_file = input_string.format(**wildcards,ext='{ext}')
            extensions, = glob_wildcards(in_file)
        
            for ext in extensions:
                src_file = in_file.format(ext=ext)
                dest_file = output_string[seq].format(**wildcards,run=run,ext=ext)
                print(f'copying {src_file} to {dest_file}')
                os.makedirs(os.path.dirname(dest_file), exist_ok=True)
                shutil.copy(src_file,dest_file)


    
    
