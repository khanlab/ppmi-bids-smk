configfile: 'config.yml'


def get_raw_zip_dir(wildcards):
    return checkpoints.download_from_csv.get(**wildcards).output[0]


def get_subject_bids_dirs(wildcards):
    """ this is an input function for a target rule, so it doesn't 
    have wildcards passed on to it, so we get them by globbing the url_csv dir"""

    prefixes, = glob_wildcards('url_csv/{prefix}.csv') 
    files = []
    for prefix in prefixes:
        ck_output = checkpoints.parse_raw_zipfiles.get(prefix=prefix).output[0]
        with open(ck_output) as f:
            subjects = json.load(f)['all_subjects']
            files.extend(expand('subj_bids/{prefix}/sub-{subject}',
                        prefix=prefix,
                        subject=subjects))
    return files



rule all_bids:
    input: 
        get_subject_bids_dirs



checkpoint download_from_csv:
    """ Run this job on CBS server node, using the web browser to get the URL csv, saved into url_csv/"""
    input: 'url_csv/{prefix}.csv'
    output: directory('raw_zips/{prefix}')
    group: 'dl'
    shell: 
        'mkdir -p {output} && for url in `cat {input}`; do wget $url --directory-prefix={output}; done '


checkpoint parse_raw_zipfiles:
    """ want to create a single json that has mapping from subject id to zipfiles"""
    input: 
        zip_dir=get_raw_zip_dir
    output: 
        'raw_zips/{prefix}_mapping.json'
    group: 'dl'
    script: 
        'scripts/get_zipfile_mapping.py'


  

rule create_subj_zip:
    """ gets dicoms for a single subject from all the zipfiles, and reorganizes them using pybids, saving them into a new zip"""
    input:
        zip_dir='raw_zips/{prefix}',
        mapping = 'raw_zips/{prefix}_mapping.json'
    output:
        zip_file='subj_zips/{prefix}/{subject}.zip'
    threads: 8
    resources: 
        mem_mb=32000,
        time=30
    group: 'subj'
    script:
        'scripts/create_subject_zip.py'





rule convert_subj_bids:
    input:
        zip_file='subj_zips/{prefix}/{subject}.zip',
        heuristic='resources/heuristic.py'
    output:
        bids_dir=directory('subj_bids/{prefix}/sub-{subject}')
    shadow:
        'minimal'
    container: config['singularity']['heudiconv']
    threads: 8
    resources: 
        mem_mb=32000,
        time=30
    group: 'subj'
    shell:
        'tmpdir=`mktemp -d` && '
        'unzip -d $tmpdir {input.zip_file} && '
        # hacky way of globbing to get sessions from the extracted zip, but avoids having to include session as a snakemake wildcard!
        'for session in $(find $tmpdir/PPMI/{wildcards.subject} -maxdepth 1 -mindepth 1 -type d| xargs -I {{}} basename {{}}); '            
        ' do '
        '  heudiconv ' 
        '  -d $tmpdir/PPMI/{{subject}}/{{session}}/*/*dcm '
        '  -s {wildcards.subject} -ss ${{session}} '
        '  --outdir $tmpdir/bids_{wildcards.subject} '
        '  --heuristic {input.heuristic} '
        '  --bids --minmeta --overwrite; '
        ' done && '
        ' mkdir -p subj_bids/{wildcards.prefix} && '
        ' cp -Rv $tmpdir/bids_{wildcards.subject}/sub-{wildcards.subject} {output.bids_dir} '
        ' && rm -rf $tmpdir '

rule get_sequences_json:
    input:
        zip_dir=get_raw_zip_dir
    output:
        'raw_zips/{prefix}_sequences.json'
    group: 'dl'
    script: 
        'scripts/get_zipfile_sequences.py'

