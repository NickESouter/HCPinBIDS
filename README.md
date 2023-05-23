# HCPinBIDS

This Python script was written by Nick Souter (N.Souter@sussex.ac.uk) in order to convert raw unprocessed data from the human connectome project (HCP) young adult dataset into brain imaging data structure (BIDS) format. 
This script is written such that it works with the file tree generated when downloading unprocessed MRI data from the ConnectomeDB portal (https://db.humanconnectome.org). For exmaple, your input directory should contain overarching folders labelled as such:

- 100307_3T_Diffusion_unproc
- 100307_3T_rfMRI_REST1_unproc
- 100307_3T_Structural_unproc
- 100307_3T_tfMRI_EMOTION_unproc
- 125525_7T_Diffusion_unproc
- 125525_7T_rfMRI_REST1_unproc
- 131722_7T_tfMRI_RETEXP_unproc

The exact folders will vary based on the participants and scans you've downloaded. It should convert all structural, functional, fieldmap, and diffusion scans, and generate JSON sidecars with metadata for each. The script will print the filename of each file being converted, while running.

## Usage

This script will expect the following arguments:

```
HCPinBIDS.py -in <input directory> -out <output directory> -method <'copy'/'link'>
```  
  
 - '-in' must be a valid directory, and should contain unprocessed HCP data in the format provided by ConnectomeDB (see above). Script will exit if an invalid path is provided.
 
 - '-out' should point to where you want your output files to be placed.
 
 - '-method' can be either 'copy' or 'link'. 'copy' will create new versions of original files in the output directory, while 'link' will create symbolic links to them with new names but without generating new files 
 (for the majority of files. Sidecar JSON files and structural magnitude fieldmaps will need to be generated as new files). In either case, files in the input directory will not be impacted.
 Note that the resuling 'link' directory will not be recognised as BIDS valid through the BIDS Validator (https://bids-standard.github.io/bids-validator), but a 'copy' version will.
 
 These arguments can be presented in any order.

## Demographic information

This script also expects a CSV file containing the participant ID, age, and gender of all participants being converted. I have included a CSV file with this information that should allow the
script to automatically pull demographics for subjects you are converting. If this file isn't present in the same directory as the script, it'll print a warning rather than crashing.

This file contains data made available by HCP under their Open Access Data Use Terms (www.humanconnectome.org/study/hcp-young-adult/document/wu-minn-hcp-consortium-open-access-data-use-terms).
As such, users must abide by these terms when using this demographic information.

## Data sources and background

Files are structured to be compatible with BIDS using the format described in v1.8.0 of the BIDS specification (https://zenodo.org/record/7263306). All available BIDS metadata for each file/scan
type was taken manually from documents provided by HCP, including Appendix 1 of the WU-Minn HCP 1200 Subjects Release: Reference Manual, and the HCP 7T Imaging Protocol Overview:

- www.humanconnectome.org/storage/app/media/documentation/s1200/HCP_S1200_Release_Reference_Manual.pdf
- www.humanconnectome.org/storage/app/media/documentation/s1200/HCP_S1200_Release_Appendix_I.pdf
- www.humanconnectome.org/hcp-protocols-ya-7t-imaging

While efforts have been made to make sure this data is as accurate and comprehensive as possible, users may consider double checking the data for any fields they are particularly concerned about, within the script.

This script was written and tested in Linux. A sample of 3T data has been succesfully run in fMRIPrep (7T data lacks the necessary T1w scans).

## Not working?

Do you download unprocessed HCP data in a way that is not compatible with this script? Do metadata files appear to be incorrectly written? Is the conversion working, but providing data that is not recognised as BIDS valid?
If so, let me know! I'll do my best to adjust the code.

It should be noted that some warnings on the BIDS validator do not reflect unsuccesful conversion, but are caused by quirks of the data provided by HCP. These include and may not be limited to:

- [Code 25] EVENTS_TSV_MISSING
- [Code 85] SUSPICIOUSLY_LONG_EVENT_DESIGN
- [Code 97] MISSING_SESSION
