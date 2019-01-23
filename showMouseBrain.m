clear all; close all; clc

load ../data/tesselated_atlas/Tesselated_Atlas.mat

faces=T;
vertices=r;
colors=cond;

figure;
patch('vertices',vertices,'faces',T,'Cdata',colors);