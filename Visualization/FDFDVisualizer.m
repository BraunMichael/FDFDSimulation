clear all; close all; clear classes; clear opts; clc;

%% Import data
addpath(genpath('E:\Michael\Stanford\Research\Data\Simulation'))
cd 'E:\Michael\Stanford\Research\Data\Simulation'
%addpath(genpath('C:\Spectre Working Folder\Reflectometry'))
%cd 'C:\Spectre Working Folder\Reflectometry'
folderpath = uigetdir;

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(folderpath)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', folderpath);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(folderpath, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
clear opts
opts.withobjsrc = true;
opts.withabs = true;
z_location = 0;

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(folderpath, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    load(fullFileName)
    figure    
    vis2d(E{Axis.x}, Axis.z, z_location, obj_array, src_array, opts)
    drawnow; % Force display to update immediately.
end


