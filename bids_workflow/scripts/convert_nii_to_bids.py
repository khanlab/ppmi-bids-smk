from snakemake.io import glob_wildcards
from snakemake.io import strip_wildcard_constraints
from snakemake.shell import shell
import os
import shutil
import re
from glob import glob
from snakebids import bids

#get subject and sessions
subject = snakemake.wildcards.subject
sessions, = glob_wildcards(f'{snakemake.input.nii_dir}/ses-{{session,[a-zA-Z0-9]+}}')

#now loop over sessions:
for session in sessions:

    #search pattern for nii.gz files in the session folder (output from dcm2niix)
    input_string_with_constraints = f'{snakemake.input.nii_dir}/ses-{session}/sub-{subject}_ses-{session}_series-{{desc,[^.]+}}.{{ext,nii.gz}}'

    #glob to get the wildcards (series description and extension)
    descs,exts = glob_wildcards(input_string_with_constraints)

    #for formatting:
    input_string = strip_wildcard_constraints(input_string_with_constraints)

    bids_dir=snakemake.params.bids_dir
    output_string = { 
                't1': bids(root=bids_dir,
                            subject=subject,
                            session=session,
                            datatype='anat',
                            run='{run}',
                            suffix='T1w.{ext}'),
                'pd': bids(root=bids_dir,
                            subject=subject,
                            session=session,
                            datatype='anat',
                            run='{run}',
                            suffix='PDw.{ext}'),
                'nm': bids(root=bids_dir,
                            subject=subject,
                            session=session,
                            datatype='anat',
                            run='{run}',
                            suffix='MTw.{ext}'),
                'flair': bids(root=bids_dir,
                            subject=subject,
                            session=session,
                            datatype='anat',
                            run='{run}',
                            suffix='FLAIR.{ext}'),
                't2': bids(root=bids_dir,
                            subject=subject,
                            session=session,
                            datatype='anat',
                            run='{run}',
                            suffix='T2w.{ext}'),
                'dwi': bids(root=bids_dir,
                            subject=subject,
                            session=session,
                            datatype='dwi',
                            run='{run}',
                            suffix='dwi.{ext}'),
                'bold': bids(root=bids_dir,
                            subject=subject,
                            session=session,
                            datatype='func',
                            task='rest',
                            run='{run}',
                            suffix='bold.{ext}'),
                }

    #initialize the mapping, indexing by session as well so we can count the runs within a session
    seq_dict= { seq:[] for seq in output_string.keys() }
        
    for desc,ext in zip(descs,exts): 

        in_file = input_string.format(desc=desc,ext=ext)
        

        descl=desc.lower() #for case-insensitive
        if ('t1' in descl) or ('mprage' in descl) or ('spgr' in descl) or ('bravo' in descl) or ('mpr' in descl and 'sag' in descl):
            seq_dict['t1'].append({'desc':desc})

        elif (('dti' in descl) or ('diff' in descl) or ('gated' in descl) or ('dwi' in descl)) and ('trace' not in descl):

            #check if bvec and bval exist:
            if os.path.exists(input_string.format(desc=desc,ext='bvec')) and os.path.exists(input_string.format(desc=desc,ext='bval')):
                seq_dict['dwi'].append({'desc': desc})

        elif 'flair' in descl:
            seq_dict['flair'].append({'desc':desc})

        elif (('rs' in descl) or ('rest' in descl)) and ('transverse' not in descl): #transverse can falsely match rs
            seq_dict['bold'].append({'desc':desc})

        elif ('nm' in descl) or ('mt' in descl):
            seq_dict['nm'].append({'desc':desc})

        elif ('pd' in descl) or ('dual' in descl):
            #check echo number if it is present:
            echo_match = re.search(r'e([0-9]+)$',desc)
            if echo_match is not None:
                echo = echo_match.group(1)
                seq_dict['t2'].append({'desc':desc})
            else:
                seq_dict['pd'].append({'desc':desc})

#        elif 't2' in descl:
#            seq_dict['t2'].append({'desc':desc})
            
#        elif 'swan' in descl:
#            seq_dict['swan'].append({'desc':desc})

        else:
            print('Cannot find matching sequence for sub-{subject}, ses-{session}, series-{desc}')


    #loop through types of sequences:
    for seq in seq_dict.keys():
        #loop through runs:
        for run,wildcards in enumerate(seq_dict[seq],start=1):

            #glob to get all the extensions for a scan:
            in_file = input_string.format(**wildcards,ext='{ext}')
            extensions, = glob_wildcards(in_file)
        
            for ext in extensions:
                src_file = in_file.format(ext=ext)
                dest_file = output_string[seq].format(**wildcards,run=run,ext=ext)
                print(f'copying {src_file} to {dest_file}')
                os.makedirs(os.path.dirname(dest_file), exist_ok=True)
                shutil.copy(src_file,dest_file)



#copy files to validate the subject - note these should only exist in the
# shadow dir, and won't conflict with the actual output files 
for bids_file in snakemake.input.bids_files:
    shutil.copy(bids_file,snakemake.params.bids_dir) 

#run the bids validator
shell('singularity run -e {snakemake.input.validator} {snakemake.params.bids_dir} --verbose > {snakemake.log}' )

