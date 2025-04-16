1.	Overview
smFISH_batchRun: A smFISH image processing tool for single-molecule RNA Detection and 3D reconstruction

    •	smFISH_batchRun is a set of custom MATLAB scripts that allow us to detect both nuclear nascent transcripts at the active transcription sites (ATS) and mature cytoplasmic mRNAs with single-molecule precision and reconstruct the tissue in 3D for further analysis. 

  •	smFISH_batchRun consists of multiple implemented functions to aid in detection and two functions that aid in visualizing the detection

  •	 smFISH_batchRun has no limit with the size of the Z-stack

  •	These codes were initially optimized for the C. elegans germline but were designed to be broadly applicable to other species and tissue types.

  •	Requirements: MATLAB (version 2019a or above)

2.	System requirement

  •	Any Windows OS or Mac OS system that can run MATLAB scripts. If using a Mac OS system, be sure to use the right input folder path in the line 50 of the function “Batch_smFISH_v2_10” before use.

  •	To run the code, the 'bioformats_package.jar' file from OpenMicroscopy is required. Make sure that the ‘bioformats_package.jar’ file is located in the following directory: ‘smFISH\Functions\bfmatlab’. The jar file can be downloaded from https://downloads.openmicroscopy.org/bio-formats/8.1.1/artifacts/bioformats_package.jar or https://www.openmicroscopy.org/bio-formats/downloads/

3.	How to use

  a.	Open the script “smFISH_batchRun” on MATLAB

  b.	Change the variable “inputPath” as desired to designate the Input folder with the images that need to be ran through the detection pipeline 

  c.	Change the variable “masterOutputPath” as desired to designate the output folder with the images that need to be ran through the detection pipeline 

  d.	Designate the thresholds greater than or equal to zero for intron and exon as desired

  e.	Specify the ChOrder as described:

    -	Put '1' for nucleus, '2' for exon, '3' for intron, '4' for protein IF, '5' for DTC channel and '0' for skipping that channel. 

  f.	Press “Run” to execute the detection pipeline in the script and generate the ‘af’ variable

  g.	Load the variable af into the MATLAB workspace

  h.	Use the af variable in the function “VisualRNAdetect” and “VisualDNAdetect” to visualize the detection of RNA and nuclei 

  i.	Repeat steps d-h as needed for obtaining the correct thresholds in each image 

This project is licensed under the terms of the MIT license
