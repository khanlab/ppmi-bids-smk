# ppmi-bids-smk

Downloading & BIDS conversion snakemake workflow for the
[PPMI](https://www.ppmi-info.org/) dataset. 

The pypmi workflow requires you to download and extract all the PPMI data
to a central location first, and it reorganizes the dicoms in-place, then 
runs heudiconv using Docker. This is not suitable for HPC systems with
inode (file number) constraints, so this workflow instead deals with 
conversion at a subject-level. Heudiconv also chokes on some dicom images 
because of incompatibility in dcmstack, so this workflow does not use
heudiconv, and instead implements the heuristic in the `convert_nii_to_bids` 
rule, and the [script](bids_workflow/scripts/convert_nii_to_bids.py) it uses.

The overall steps of the workflow are:
 - Downloading the zip files
 - Parsing the list of files in the zips without extracting, to
  obtain the list of subjects, and which zipfiles are associated with 
  each subject
 - Extracting the dicom files for each subject into a /tmp folder,
   reorganizing them, and saving as an individual zip
 - Creating the bids subject dir by extracting each subject zip 
  into a /tmp folder, converting to nifti with dcm2niix, and copied 
  into bids
 
The "heuristic" for BIDS conversion is now embedded in this 


The workflow is split into two, one to download and re-package 
the zipfiles, and another to convert those to bids. Combined together,
this is an end-to-end fully-automated workflow.

## Instructions

 1. Install snakemake and pypmi 
 2. Use a web browser (see note 1) to download the URL csv file into 
   the `url_csv` folder
 3. Run a dry-run of the `zips_workflow` to download (dry-run, print shell cmd, reason)::
    cd zips_workflow
    snakemake -npr --until download_from_csv
 4. Actually run it::
    snakemake --cores all --until download_from_csv
 5. Dry-run the rest of the `zips_workflow` (see note 2)::
    snakemake -npr
 6. Run the entire `zips_workflow` to create the `subj_zips` folder::
    snakemake --cores all
 7. Now run the `bids_workflow` to create the `subj_bids` folder::
    cd ../bids_workflow
    snakemake -npr
    snakemake --cores all



 Note 1: The IP of the system used for steps 2-4 (creating/downloading the URL csv file, 
  and for downloading the zip files) must be identical. 
  Since you want to run the workflow on a HPC that does not have a browser, one way to 
  do this is to use `sshfs` to mount the remote server on the system used for downloading.
  Or you can just copy the data afterwards. 

 Note 2: Make sure this step does not try to run the `download_from_csv` rule again,
   which would end up clearing your download directory. You may need to use the following 
   command to clean-up the metadata to avoid this
   
    snakemake --cleanup-metadata raw_zips/*

