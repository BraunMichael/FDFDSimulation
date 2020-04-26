close all;
clear all;
clc;
%Select folder containing all the files
folderpath=uigetdir();
%Set this folder as current directory
cd(folderpath);


% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(folderpath)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', folderpath);
    uiwait(warndlg(errorMessage));
    return;
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(folderpath, '*.mat');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    Simfilename = fullfile(folderpath, baseFileName);
    [~, Simfilename_only, ~]=fileparts(Simfilename);
    
    % %% Import data from text file.
    % addpath(genpath('E:\Michael\Stanford\Research\Data\Simulation'))
    % cd 'E:\Michael\Stanford\Research\Data\Simulation'
    % %addpath(genpath('C:\Spectre Working Folder\Reflectometry'))
    % %cd 'C:\Spectre Working Folder\Reflectometry'
    % [Simfilename, folderpath] = uigetfile('*.mat');
    % cd(folderpath)
    %[~, Simfilename_only, ~]=fileparts(Simfilename);
    load(Simfilename)
    fprintf('Working on %s\n', Simfilename_only);
    Fields = {E, H};
    length = obj_array(2).shape.L(2);
    diameter = obj_array(2).shape.L(1);
    i=1;
    for field = Fields
        for component = Axis.elems
            fprintf('Working on %s\n', field{1}{component}.name)
            primx = field{1}{component}.grid3d.lall{1,1};
            primy = field{1}{component}.grid3d.lall{2,1};
            primz = field{1}{component}.grid3d.lall{3,1};
            dualx = field{1}{component}.grid3d.lall{1,2};
            dualy = field{1}{component}.grid3d.lall{2,2};
            dualz = field{1}{component}.grid3d.lall{3,2};
            
            interpx = sort(horzcat(primx,dualx));
            interpy = sort(horzcat(primy,dualy));
            %interpz = sort(horzcat(primz,dualz));
            %Only calculate at z=0
            interpz = 0;
            
            [Xq, Yq, Zq] = meshgrid(interpx, interpy, interpz);
            switch field{1}{component}.gt_array(Axis.x)
                case 'prim'
                    tempX = primx;
                case 'dual'
                    tempX = dualx;
            end
            switch field{1}{component}.gt_array(Axis.y)
                case 'prim'
                    tempY = primy;
                case 'dual'
                    tempY = dualy;
            end
            switch field{1}{component}.gt_array(Axis.z)
                case 'prim'
                    tempZ = primz;
                case 'dual'
                    tempZ = dualz;
            end
            [X, Y, Z] = meshgrid(tempX, tempY, tempZ);
            
            %Use permute to switch from Maxwell_FDFD to Matlab standard of Y,X,Z array
            %indices
            field_values = permute(field{1}{component}.array, [2 1 3]);
            
            interpfield_values = griddata(X, Y, Z, field_values, Xq, Yq, Zq, 'natural');
            Compiled_fields{i}{component} = interpfield_values;
        end
        i = i + 1;
    end
    E_compiled = cat(3, Compiled_fields{1}{Axis.x}, Compiled_fields{1}{Axis.y}, Compiled_fields{1}{Axis.z});
    H_compiled = cat(3, Compiled_fields{2}{Axis.x}, Compiled_fields{2}{Axis.y}, Compiled_fields{2}{Axis.z});
    
    for i=1:size(E_compiled,2)
        for j=1:size(E_compiled,1)
            E2(j, i) = sqrt(E_compiled(j, i, 1)*conj(E_compiled(j, i, 1)) + E_compiled(j, i, 2)*conj(E_compiled(j, i, 2)) + E_compiled(j, i, 3)*conj(E_compiled(j, i, 3)));
            H2(j, i) = sqrt(H_compiled(j, i, 1)*conj(H_compiled(j, i, 1)) + H_compiled(j, i, 2)*conj(H_compiled(j, i, 2)) + H_compiled(j, i, 3)*conj(H_compiled(j, i, 3)));
        end
    end
    
    S = cross(E_compiled,conj(H_compiled));
    
    for i=1:size(S,2)
        for j=1:size(S,1)
            S_mag(j, i) = sqrt(S(j, i, 1)*conj(S(j, i, 1)) + S(j, i, 2)*conj(S(j, i, 2)) + S(j, i, 3)*conj(S(j, i, 3)));
        end
    end
    
    save(Simfilename, 'E', 'H', 'obj_array', 'src_array', 'J', 'E_compiled', 'H_compiled', 'E2', 'H2', 'S', 'S_mag', 'Xq', 'Yq', 'Zq')
end

fprintf('Done')
