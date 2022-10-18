from snakemake.shell import shell
from glob import glob
import tempfile
import json
import xmltodict


subject=snakemake.wildcards.subject

with tempfile.TemporaryDirectory() as tmpdirname:
 
    #unzip to tmpdir
    shell('unzip -d {tmpdirname} {snakemake.input.zip_file}')

    #loop over image directories - the bottom-most folder is the unique identifier
    # we will use to grab the metadata info on the subject 
        

    for session_dir in glob(f'{tmpdirname}/PPMI/{subject}/*/*'):

        study_id = session_dir.split('/')[-1] #this is unique for each scan

        #pull out the corresponding metadata xml from the metadata zip
        shell('unzip -d {tmpdirname}/metadata {snakemake.input.metadata_zip} PPMI/PPMI_{subject}_*_{study_id}.xml')

        #now read that xml file in, using xmltodict to make it easier to use
        xml_file = glob(f'{tmpdirname}/metadata/PPMI/PPMI_{subject}_*_{study_id}.xml')[0]
        with open(xml_file,'r') as f:
            metadata_dict = xmltodict.parse(f.read())


        metadata_dict = metadata_dict['idaxs']['project'] #leave out top-levels

        #replace @term, #text with simple dict instead for protocolTerm
        for protocol_term in metadata_dict['subject']['study']['imagingProtocol']['protocolTerm']['protocol']:
            metadata_dict['subject']['study']['imagingProtocol'][protocol_term['@term']] = protocol_term.get('#text','')
        del metadata_dict['subject']['study']['imagingProtocol']['protocolTerm']



        visit = metadata_dict['subject']['visit']['visitIdentifier']
        session = visit.replace(' ','').replace('_','').replace('-','') #remove whitespace, _ and - for bids

        series_desc = metadata_dict['subject']['study']['imagingProtocol']['description']
        series_desc = series_desc.replace(' ','_') #replace space with underscores

        #could use bids naming at this point already... but we will save it for another rule
        shell("mkdir -p {snakemake.output.nii_dir}/ses-{session}")


        print(f'converting dicoms from {session_dir} as sub-{subject}_ses-{visit}_series-{series_desc}')
        shell("singularity exec -e {snakemake.params.container_uri} dcm2niix -o {snakemake.output.nii_dir}/ses-{session} -f 'sub-{subject}_ses-{session}_series-{series_desc}' {session_dir}")

        #add the ppmi metadata to the json file(s)
        json_files = glob(f'{snakemake.output.nii_dir}/ses-{session}/sub-{subject}_ses-{session}_series-{series_desc}*json')
        for json_file in json_files:

            print(f'updating {json_file}')
            with open(json_file,'r') as f:
                try:
                    json_dict = json.load(f)
                except:
                    print(f'WARNING: dcm2niix produced empty json for {json_file}')
                    json_dict = {'convert_zip_to_nii':'dcm2niix produced empty json'}

            json_dict['SequenceVariant'] = ""

            json_dict['ppmi'] = metadata_dict
            with open(json_file,'w',encoding='utf-8') as f:
                json.dump(json_dict, f, indent=4)





