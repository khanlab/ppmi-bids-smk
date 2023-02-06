import json
from snakemake.shell import shell
from pypmi import bids
import tempfile


with open(snakemake.input.mapping) as f:
    mappings = json.load(f)

zipfiles = mappings['subjects_to_zipfiles'][snakemake.wildcards.subject]
print(zipfiles)

with tempfile.TemporaryDirectory() as tmpdir:

    #unzip the subject files from each zipfile
    for f in zipfiles:
        shell(
            "unzip -o {f}  PPMI/{snakemake.wildcards.subject}/* -d {tmpdir} "
            )

    #run pypmi clean-up on the files
    if snakemake.params.pypmi_prepare:
        subjects, coerced = bids._prepare_subject(f'{tmpdir}/PPMI/{snakemake.wildcards.subject}',confirm_uids=False)


        print(subjects)
        print(coerced)

    #zip those and copy back
    shell(
        "pushd {tmpdir} && "
        "zip -0 {snakemake.wildcards.subject}.zip -r PPMI/{snakemake.wildcards.subject} && "
        "popd && "
        "cp -v {tmpdir}/{snakemake.wildcards.subject}.zip {snakemake.output.zip_file}"
        )
        

