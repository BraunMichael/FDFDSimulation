clear all; close all; clear classes; clear opts; clc;

set(0,'defaultaxesfontname','cmu serif')
set(0,'DefaultAxesFontSize',14)

folderpath=uigetdir();
addpath(genpath(folderpath))
cd(folderpath)

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(folderpath, '*.mat');
theFiles = dir(filePattern);
[tmp, ind] = sort({theFiles.name});
theFiles = theFiles(ind);
integ_series = zeros(1, length(theFiles));
length_series = zeros(1, length(theFiles));
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    Simfilename = fullfile(folderpath, baseFileName);
    [~, intermediate, ~] = fileparts(Simfilename);
    [~, Simfilename_only, ~] = fileparts(intermediate);
    load(Simfilename_only)
    fprintf('Working on %s\n', Simfilename_only);
    
    u_int = u(:,1)';
    v_int = v(1,:);
    
    integ_series(k) = trapz(v_int, trapz(u_int,Fh, 1));
    %length_series(k) = 1340 - 100 - L(k);
    if strfind(Simfilename_only, 'L_0_')
        length = 0;
    end
    length_series(k) = length;
end
[length_series, length_indecies] = sort(length_series);
integ_series = integ_series(length_indecies);

disp('Done')
figure();
plot(length_series,integ_series,'ko-')
title('s=300 Source - Substrate Distance Background Emission Subtracted','Interpreter','latex')
%title('s=300 vs 400 - Background Emission Subtracted','Interpreter','latex')
xlabel('$$L$$ (nm)','interpreter','latex')
ylabel('Intensity (arb units)','interpreter','latex')
axis square
%legend('7 wires/um2', '11 wires/um2')