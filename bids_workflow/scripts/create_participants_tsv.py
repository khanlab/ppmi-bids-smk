from glob import glob
import json
import pandas as pd

df_list=[]
for i,subj_dir in enumerate(sorted(snakemake.input.subj_dirs)):

    print(subj_dir)
    #glob to get the first json in the data folders
    bids_jsons = sorted(glob(f'{subj_dir}/ses-*/*/*.json'))

    for bids_json in bids_jsons:


        with open(bids_json) as f:
            json_metadata = json.load(f)

        if 'ppmi' not in json_metadata.keys():
            continue
        
            
        metadata = json_metadata['ppmi']
        subj_metadata = {'participant_id': f"sub-{metadata['subject']['subjectIdentifier']}",
                        'group': metadata['subject']['researchGroup'],
                        'age': metadata['subject']['study']['subjectAge'],
                        'weightKg': metadata['subject']['study']['weightKg'],
                        'sex': metadata['subject']['subjectSex'],
                        'site': metadata['siteKey'],
                        'visit': metadata['subject']['visit']['visitIdentifier'],
                        'date_acquired': metadata['subject']['study']['series']['dateAcquired'],
                        'scanner_mfg': metadata['subject']['study']['imagingProtocol']['Manufacturer'],
                        'scanner_make': metadata['subject']['study']['imagingProtocol']['Mfg Model'],
                        'field_strength': metadata['subject']['study']['imagingProtocol']['Field Strength']}

        df_subj=pd.DataFrame(subj_metadata,index=[i])

        df_list.append(df_subj)
        break


df = pd.concat(df_list)
print(df)
df.to_csv(snakemake.output.tsv,sep='\t',index=False)
