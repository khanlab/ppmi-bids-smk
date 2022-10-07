# ppmi-bids-smk

Downloading & BIDS conversion snakemake workflow. Makes use of 
[pypmi](https://github.com/rmarkello/pypmi) for dicom reorganization 
and for the PPMI BIDS conversion heuristic file. 

The pypmi workflow requires you to download and extract all the PPMI data
to a central location first, and it reorganizes the dicoms in-place, then 
runs heudiconv using Docker. This is not suitable for HPC systems with
inode (file number) constraints, so this workflow instead deals with 
conversion at a subject-level. 

It does this by:
 - Downloading the zip files
 - Parsing the list of files in the zips without extracting, to
  obtain the list of subjects, and which zipfiles are associated with 
  each subject
 - Extracting the dicom files for each subject into a /tmp folder,
   reorganizing them, and saving as an individual zip
 - Creating the bids subject dir by extracting each subject zip 
  into a /tmp folder, then converting to BIDS with heudiconv
   
## Instructions

 1. Install snakemake and pypmi 
 2. Use a web browser (see note 1) to download the URL csv file into 
   the `url_csv` folder
 3. Run a dry-run of the workflow to download (prints what it will do)::
    snakemake --dry-run --until download_from_csv
 4. Actually run it::
    snakemake --cores all --until download_from_csv
 5. Run the rest of the workflow to convert::
    snakemake --cores all


 Note 1: The IP of the system used for steps 2-4 (creating/downloading the URL csv file, 
  and for downloading the zip files) must be identical. 
  Since you want to run the workflow on a HPC that does not have a browser, one way to 
  do this is to use `sshfs` to mount the remote server on the system used for downloading.
  Or you can just copy the data afterwards. 
