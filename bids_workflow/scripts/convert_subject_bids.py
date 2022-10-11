from snakemake.shell import shell
from pypmi import bids
import tempfile

with tempfile.TemporaryDirectory() as tmpdirname:
    

    #unzip into the shadow dir
    shell("unzip -d {tmpdirname}  {snakemake.input.zip_file}")


    subjects, coerced = bids._prepare_subject(raw_dir=f'{tmpdirname}/PPMI/{snakemake.wildcards.subject}',ignore_bad=True,coerce_study_uids=False)
    

    shell("zip -r {tmpdirname}  {snakemake.input.zip_file}")
                     
