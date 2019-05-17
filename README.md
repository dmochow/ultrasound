# ultrasound
This is the code repository for the Dmochowski Lab's ultrasonic stimulation project.  

You must have the data on your machine with the following folder structure:
data/r0xx/activeYY/filename.mat

where xx is the rat number (e.g. 012) and YY is either '10' or '40'
filename must indicate whether this is the baseline or post-stimulation recording

runHilbert_debug.m preprocesses the data of each rat/condition (separately) and saves the data into a mat file
analyzeHilbert.m is a group analysis of the preprocessed data

(c) 2018- Jacek "Loquacious D" Dmochowski
