% Script to analyze data of particle length measurements
clear all; close all;

ss = dir('*.mat');

% Load all particle file data
for k = 1:length(ss)
    Q{k} = load(ss(k).name);
    Q{k}.filename = ss(k).name; 
end

% Join all data into single dataset.
particlelenum = [];
particlenum = [];
particlefile = [];

for kfile = 1:length(Q)
    for kpart = 1:length( Q{kfile}.particles )
        particlelenum(end+1) = Q{kfile}.particles(kpart).totallen_um;
        particlenum(end+1) = kpart;
        particlefile{end+1} = Q{kfile}.filename;
    end
end
figure;
hist(particlelenum,sqrt(length(particlenum)))