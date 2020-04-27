%% Field Visualization
close all;
clear all;
%% Import data from mat file.
addpath(genpath('E:\Michael\Stanford\Research\Data\Simulation'))
cd 'E:\Michael\Stanford\Research\Data\Simulation'
%addpath(genpath('C:\Spectre Working Folder\Reflectometry'))
%cd 'C:\Spectre Working Folder\Reflectometry'
[Simfilename, folderpath] = uigetfile('*.mat');
cd(folderpath)
[~, Simfilename_only, ~]=fileparts(Simfilename);
clc;
load(Simfilename)
fprintf('Working on %s\n', Simfilename_only);
Fields = {E, H};
Compiled_Fields = {E_compiled, H_compiled};
%labels = {'E_x', 'E_y', 'E_z', 'H_x', 'H_y', 'H_z'}
length = obj_array(2).shape.L(2);
diameter = obj_array(2).shape.L(1);
i = 1;
j = 1;
for field = Fields
    for component = Axis.elems
        fprintf('Working on %s', field{1}{component}.name)
        figure
        axhand = gca;
        fighand = pcolor(Xq, Yq, abs(interpfield_values));
        daspect(axhand, [1 1 1])
        set(fighand, 'EdgeColor', 'none')
        set(fighand, 'FaceColor', 'interp')
        title(sprintf('%s',field{1}{component}.name),'Interpreter','latex')
        xlabel('Location (nm)','Interpreter','latex')
        ylabel('Location (nm)','Interpreter','latex')
        colormap(axhand, jet(4096))
        colorbar
        drawnow
        j = j + 1;
    end
    i = i + 1;
end

figure
axhand = gca;
fighand = pcolor(Xq, Yq, abs(E2));
daspect(axhand, [1 1 1])
set(fighand, 'EdgeColor', 'none')
set(fighand, 'FaceColor', 'interp')
title('E2','Interpreter','latex')
xlabel('Location (nm)','Interpreter','latex')
ylabel('Location (nm)','Interpreter','latex')
colormap(axhand, jet(4096))
colorbar

S_labels = {'S_x', 'S_y', 'S_z'};
for i = 1:max(size(Slabels))
figure
axhand = gca;
fighand = pcolor(Xq, Yq, abs(S(:,:,i)));
daspect(axhand, [1 1 1])
set(fighand, 'EdgeColor', 'none')
set(fighand, 'FaceColor', 'interp')
title(S_labels{i},'Interpreter','latex')
xlabel('Location (nm)','Interpreter','latex')
ylabel('Location (nm)','Interpreter','latex')
colormap(axhand, jet(4096))
colorbar
end

figure
axhand = gca;
fighand = pcolor(Xq, Yq, abs(S_mag));
daspect(axhand, [1 1 1])
set(fighand, 'EdgeColor', 'none')
set(fighand, 'FaceColor', 'interp')
title('S_mag','Interpreter','latex')
xlabel('Location (nm)','Interpreter','latex')
ylabel('Location (nm)','Interpreter','latex')
colormap(axhand, jet(4096))
colorbar