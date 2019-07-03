# ultrasound
This is the code repository for the Dmochowski Lab's ultrasonic stimulation project.  

You must have the data on your machine with the following folder structure:
data/r0xx/activeYY/filename.mat

where xx is the rat number (e.g. 012) and YY is either '10' or '40'
filename must indicate whether this is the baseline or post-stimulation recording


Procedures:

(1) callReadIntanFolder.m: to convert RHD files into mat files (still raw data)
(2) runHilbert_debug.m: only up to line 142 (feel free to rip this code and make new function for preproc)
(3) analyzeHilbert.m, analyzeHilbertWithStim.m is my actual analysis



runHilbert_debug.m preprocesses the data of each rat/condition (separately) and saves the data into a mat file
analyzeHilbert.m is a group analysis of the preprocessed data

(c) 2018- Jacek "Loquacious D" Dmochowski
