# HCPinBIDS

This Python script was written by Nick Souter (N.Souter@sussex.ac.uk) in order to convert raw unprocessed data from the human connectome project (HCP) young adult dataset into brain imaging data structure (BIDS) format. 
This script is written such that it works with the file tree generated when downloading unprocessed MRI data from the ConnectomeDB portal (https://db.humanconnectome.org).
It should convert all structural, functional, fieldmap, and diffusion scans.

## Running the script

Upon running the script, you will have the option to 'copy' files into new BIDS files, or to symbolic 'link' them with new names but without creating new files.
Note that the resuling 'link' directory will not be recognised as BIDS valid through the BIDS Validator (https://bids-standard.github.io/bids-validator), but the 'copy' version will.

As it stands, the output directory is generated within the current working directory holding this script. The input_dir variable should be updated as needed within the script.

## Demographic information

This script also expects a CSV file containing the participant ID, age, and gender of all participants being converted. I have included a CSV file with this information that should allow the
script to automatically pull demographics for subjects you are converting. If this file isn't present, it'll print a warning rather than crashing.

This file contains data made available by HCP under their Open Access Data Use Terms (www.humanconnectome.org/study/hcp-young-adult/document/wu-minn-hcp-consortium-open-access-data-use-terms).
As such, users must abide by these terms when using this demographic information.

## Data sources and background

Files are structured to be compatible with BIDS using the format described in v1.8.0 of the BIDS specification (https://zenodo.org/record/7263306#.ZGdUa3bMKUk). All available metadata for each file/scan
type was taken manually from documents provided by HCP, including Appendix 1 of the WU-Minn HCP 1200 Subjects Release: Reference Manual, and the HCP 7T Imaging Protocl Overview:

www.humanconnectome.org/storage/app/media/documentation/s1200/HCP_S1200_Release_Reference_Manual.pdf
www.humanconnectome.org/storage/app/media/documentation/s1200/HCP_S1200_Release_Appendix_I.pdf
www.humanconnectome.org/hcp-protocols-ya-7t-imaging

While efforts have been made to make sure this data is as accurate and comprehensive as possible, users may consider double checking the data for any fields they are particularly concerned about, within the script.

This script was written and tested in Linux.

## Not working?

Do you download unprocessed HCP data in a way that is not compatible with this script? Do metadata files appear to be incorrectly written? Is the conversion working, but providing data that is not recognised as BIDS valid?
If so, let me know! I'll do my best to adjust the code.

It should be noted that some warnings on the BIDS validator do not reflect unsuccesful conversion, but are caused by quirks of the data provided by HCP. These include and may not be limited to:

[Code 25] EVENTS_TSV_MISSING

[Code 85] SUSPICIOUSLY_LONG_EVENT_DESIGN

[Code 97] MISSING_SESSION
