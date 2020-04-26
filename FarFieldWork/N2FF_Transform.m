%% Near field to far field transformation
close all; clear all; clc;
tic
remote = false;
limited_range = true; %limits angular range, increases point resolution

if limited_range
    center_theta = pi/6;
    center_phi = pi;
    anglerange_deg = 6.6; %Sets range of angles used if limited_range is true, corresponds to actual angles
    points_per_deg = 10; %Sets number of points per degree
else
    n_theta = 360;
    n_phi = 360;
end

if remote
    %For use remotely
    folderpath = '$HOME/mbraun7';
    addpath(genpath(folderpath))
    folderpath = '$HOME/Maxwell_FDFD';
    addpath(genpath(folderpath))
    folderpath = '$HOME/NtoFField_packageV8/';
    addpath(genpath(folderpath))
    %Set folder containing all the files
    folderpath = strcat(folderpath, 'FieldData');
    % Check to make sure that folder actually exists.  Warn user if it doesn't.
    if ~isdir(folderpath)
        fprintf('Error: The following folder does not exist:\n%s', folderpath);
        return;
    end
else %at home
    folderpath = uigetdir();
    cd(folderpath)
end


%% Setup box
% specific parameters (user does not change them)
option_cal_champ=2;option_i=1;

% physical parameters
ld=633e-9; % wavelength in m
k0=2*pi/ld;
n_strates=[1, 5.483850];
z_strates=[0]; %Might want to read in substrate thickness and enter as z_strates
grid_spacing = 20; %point spacing in nm
distance_above = 20; %nm above source to use in integration
domainsize = 10; %in um, x-y dimensions of box to integrate


%% Field Prep

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(folderpath, '*.E.h5');
theFiles = dir(filePattern);
[tmp, ind] = natsortfiles({theFiles.name});
theFiles = theFiles(ind);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    Simfilename = fullfile(folderpath, baseFileName);
    [~, intermediate, ~] = fileparts(Simfilename);
    [~, Simfilename_only, ~] = fileparts(intermediate);
    load(Simfilename_only)
    fprintf('Working on %s\n', Simfilename_only);
    [E, H] = read_output(Simfilename_only);
    save('saved_field_matrix.mat', 'E', 'H', 'obj_array', 'src_array', '-v7.3')
    primy_temp = E{Axis.x}.grid3d.lall{2,1};
    min_y = min(min(primy_temp)); max_y = max(max(primy_temp));
    
    % Distance above source in nm
    height = src_array.intercept + distance_above;
    
    % definition of the box in COMSOL coordinates (the third coordinate, z, is perpendicular to the layers)
    Centre=[0, 0, mean([min_y height])] * 1e-9; % in meters
    L=[domainsize * 1e-6, domainsize * 1e-6, abs(height - min_y) * 1e-9]; %in meters
    N_gauss = round((L .* 1e9)./ grid_spacing) + 1;
    %% Calculate init
    init = retop(@FieldInterpolation_onbox_v3,k0,n_strates,z_strates,Centre,L,N_gauss,struct('option_cal_champ',option_cal_champ,'option_i',option_i));
    return
    %% Compute the plane wave decomposition in free space
    if limited_range
        theta = linspace(center_theta - (anglerange_deg*0.5*pi/180), center_theta +(anglerange_deg*0.5*pi/180), points_per_deg*anglerange_deg + 1);
        phi = linspace(center_phi - (anglerange_deg*0.5*pi/180), center_phi +(anglerange_deg*0.5*pi/180), points_per_deg*anglerange_deg + 1);
    else
        theta = linspace(0,pi/2,n_theta);
        phi = linspace(0,2*pi,n_phi);
    end
    % if one fixes phi, one may get the polar plot in the different planes XY, XZ, ZY (2D diagram pattern)
    [Teta,Phi]=ndgrid(theta,phi);
    u=sin(Teta).*cos(Phi);
    v=sin(Teta).*sin(Phi);
    w=cos(Teta);
    
    uh=k0*n_strates(1)*u;
    vh=k0*n_strates(1)*v;
    angles_h=retop(init,uh,vh,1,struct('test',1));% calculate top angles
    Fh=reshape(angles_h.F,size(uh));
    u_int = u(:,1)';
    v_int = v(1,:);
    integ = trapz(v_int, trapz(u_int,Fh, 1));
    length = obj_array(2).shape.L(2);
    
    savefilename = strcat(Simfilename_only, '_farfield.mat');
    save(savefilename, 'angles_h', 'init', 'Fh', 'uh', 'vh', 'u', 'v', 'w', 'integ', 'length')
    
    if remote
    else
        figure;
        surf(u.*Fh,v.*Fh,w.*Fh,'facecolor',[.5,.5,.5],'LineStyle','none');
        axis tight;axis equal;
        xlabel('x');ylabel('y');zlabel('z');
        view([-20,40]);
        lighting gouraud;light('position',[-1,0,1],'color','r');light('position',[-1,-5,1],'color','w');light('position',[-1,-5,1],'color','y');
        title('radiation diagram');
        %ret_normalise_dop(1.e13);
    end
    
end
toc
%


