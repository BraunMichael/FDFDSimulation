close all; clear all; clc;
[thefile, folderpath] = uigetfile('*.E.h5');
addpath(genpath(folderpath))
cd(folderpath)
Simfilename = fullfile(folderpath, thefile);
[~, intermediate, ~] = fileparts(Simfilename);
[~, Simfilename_only, ~] = fileparts(intermediate);
load(Simfilename_only)
fprintf('Working on %s\n', Simfilename_only);
[E_background, H_background] = read_output(Simfilename_only);

background_domain_height = obj_array.shape.lprim{2}(2);
source_height_background = src_array.intercept;
obj_array_background = obj_array;
src_array_background = src_array;
save('background_fields_400_maxL_2000_200nmtop_2x2.mat', 'E_background', 'H_background', 'obj_array_background', 'src_array_background', 'background_domain_height', 'source_height_background');

disp('Done')