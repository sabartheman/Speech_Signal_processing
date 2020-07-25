%% Skylar Tamke - Final Project
% The idea of this project is to take a bunch of different TIMIT speakers
% saying the same phrase and push their formants through a clustering
% method, resulting in means for each pheonem and synthesize the wave to
% to see how the average speaker of the samples would sound
clc
clear
close all

wavfiles = dir("*.WAV")
phnfiles = dir("*.PHN")

filename1 = "SX291.WAV";
phnfile1   = "SX291.PHN";
filename2 = "SX291.WAV";
phnfile2   = "SX291.PHN";
filename3 = "SX291.WAV";
phnfile3   = "SX291.PHN";
filename4 = "SX291.WAV";
phnfile4   = "SX291.PHN";
filename5 = "SX291.WAV";
phnfile5   = "SX291.PHN";
filename6 = "SX291.WAV";
phnfile6   = "SX291.PHN";

%% Find formants and parse into new datasets for clustering
% The formant method we working on for project 2 and parse resulting
% formants into datasets.

numSets = 6
dataset = {numSets}
charset = {numSets}
[dataset{:,:,1},charset{:,1}] = parsewaveform(filename1,phnfile1);
[dataset{:,:,2},charset{:,2}] = parsewaveform(filename2,phnfile2);
[dataset{:,:,3},charset{:,3}] = parsewaveform(filename3,phnfile3);
[dataset{:,:,4},charset{:,4}] = parsewaveform(filename4,phnfile4);
[dataset{:,:,5},charset{:,5}] = parsewaveform(filename5,phnfile5);
[dataset{:,:,6},charset{:,6}] = parsewaveform(filename6,phnfile6);

data  = []
chars = []
for i = 1:numSets
   data = vertcat(data,dataset{:,:,i})
%    chars = vertcat(chars,charset{:,i})
end
chars =  {'h#';'hh';'eh';'l';'pcl';'p';'gcl';'g';'r';'ey';'gcl';'t';'ix';'pcl';'p';'ih';'kcl';'k';'ix';'pcl';'p';'eh';'kcl';'k';'ix';'v';'pcl';'p';'ix';'tcl';'t';'ey';'dx';'ow';'z';'h#';'h#';'hh';'eh';'l';'pcl';'p';'gcl';'g';'r';'ey';'gcl';'t';'ix';'pcl';'p';'ih';'kcl';'k';'ix';'pcl';'p';'eh';'kcl';'k';'ix';'v';'pcl';'p';'ix';'tcl';'t';'ey';'dx';'ow';'z';'h#';'h#';'hh';'eh';'l';'pcl';'p';'gcl';'g';'r';'ey';'gcl';'t';'ix';'pcl';'p';'ih';'kcl';'k';'ix';'pcl';'p';'eh';'kcl';'k';'ix';'v';'pcl';'p';'ix';'tcl';'t';'ey';'dx'; ...
'ow';'z';'h#';'h#';'hh';'eh';'l';'pcl';'p';'gcl';'g';'r';'ey';'gcl';'t';'ix';'pcl';'p';'ih';'kcl';'k';'ix';'pcl';'p';'eh';'kcl';'k';'ix';'v';'pcl';'p';'ix';'tcl';'t';'ey';'dx';'ow';'z';'h#';'h#';'hh';'eh';'l';'pcl';'p';'gcl';'g';'r';'ey';'gcl';'t';'ix';'pcl';'p';'ih';'kcl';'k';'ix';'pcl';'p';'eh';'kcl';'k';'ix';'v';'pcl';'p';'ix';'tcl';'t';'ey';'dx';'ow';'z';'h#';'h#';'hh';'eh';'l';'pcl';'p';'gcl';'g';'r';'ey';'gcl';'t';'ix';'pcl';'p';'ih';'kcl';'k';'ix';'pcl';'p';'eh';'kcl';'k';'ix';'v';'pcl';'p';'ix';'tcl';'t';'ey';'dx';'ow';'z';'h#'}

%creating color coding
for i = 1:3
   a(i,:) = linspace(.2,1,12);
end

dataLength = length(data);
figure(1)
for i = 1:dataLength
    scatter(data(i,1),data(i,2))
    hold on
end
