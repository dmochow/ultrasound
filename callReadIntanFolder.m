clear all; close all; clc;  
rStr='r024';
read_Intan_RHD2000_files(['../data/' rStr '/active40']);
read_Intan_RHD2000_files(['../data/' rStr '/sham40']);
read_Intan_RHD2000_files(['../data/' rStr '/active10']);
read_Intan_RHD2000_files(['../data/' rStr '/sham10']);
% delete rhd files
% call runHilbert_debug