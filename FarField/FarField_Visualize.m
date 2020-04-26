clear all; close all; clear classes; clear opts; clc;

set(0,'defaultaxesfontname','cmu serif')
set(0,'DefaultAxesFontSize',14)
%% Simulation Setup Parameters

endstring = '_L_source_top_distance_farfield.mat';
%endstring = '_farfield.mat';

radial_growth = false;
Dim_3D = true; %3D vs 2D?

%% Simulation Domain Parameters
grid_size = 10;  % grid cell size in nm
UniformNWCells = false; % True makes xyz SIZE of grid cells the same in the nanowire, False makes xyz NUMBER of grid cells the same in the nanowires

%% Simulation Parameters

substrate_thickness = 100; %633/2; %Wavelength/2 %Was 100 for most
Laser_angle = 60; %Angle in degrees defined as positive from substrate, 90degrees is perpendicular to substrate
Laser_pol = 0; %For TFSF: Angle defined as positive from the x-y plane, coming out of the page
%For normally, defined as positive from the +z axis towards the +x axis. +z is 0

num_lengths = 51;
%num_spacings = 5;
%num_diameters = 31;
%L = linspace(0, 1000, num_lengths);
L = linspace(50, 600, 56);

%L = 1000;
%S = linspace(300, 1000, num_spacings);
%S = [300 400 600 1000];
S = 300;
%D = linspace(40, 100, num_diameters);
D = 40;
%max_NWlength = max(L); % in nm
max_NWlength = 1000;

if UniformNWCells
    typestring = 'nonuniform_';
else
    typestring = 'nonuniform_xy_';
end
if radial_growth
    typestring = strcat(typestring, 'radial_');
    if max(size(D)) == 1
        disp('radial_growth is true but you only have 1 diameter! You probably meant to change one of them')
        return
    end
end

folderpath=uigetdir();
addpath(genpath(folderpath))
cd(folderpath);
integ_series = zeros(size(L));
length_series = zeros(size(L));

for i=1:max(size(S))
    if radial_growth
        % Needs to be updated to work
        
        %         if Dim_3D
        %             folderstring = sprintf('3D_%ss_%.0f_maxL_%.0f_res_%.0f%s', typestring, S(i), max_NWlength, grid_size, endstring);
        %         else
        %             folderstring = sprintf('2D_%ss_%.0f_maxL_%.0f_res_%.0f%s', typestring, S(i), max_NWlength, grid_size, endstring);
        %         end
        %         fullfolderstring = strcat(pwd,'\',folderstring);
        %         if 7~=exist(fullfolderstring,'dir')
        %             mkdir(fullfolderstring)
        %         end
        %         addpath(genpath(fullfolderstring))
        %         cd (fullfolderstring)
        %         runfilesstring = 'RunFiles';
        %         fullrunfilesstring = strcat(pwd,'\',runfilesstring);
        %         chmodfilename = strcat(fullrunfilesstring,'\','chmodfile.sh');
        %         if 7~=exist(fullrunfilesstring,'dir')
        %             mkdir(fullrunfilesstring)
        %         end
        %         addpath(genpath(fullrunfilesstring))
        %         cd(fullrunfilesstring)
        %         fid = fopen(chmodfilename,'w');
        %         fprintf(fid, '#! /bin/bash\ncd $SCRATCH/%s/RunFiles\nchmod -R +x *\ngzip -d -r $SCRATCH/%s/\n', folderstring, folderstring);
        %         fclose(fid);true;
    end
    for j=1:max(size(D))
        for m=1:max(size(Laser_angle))
            for n=1:max(size(Laser_pol))
                for p=1:max(size(substrate_thickness))
                    for k=1:max(size(L))
                        runningfilenamebase = sprintf('GeNW_%sL_%0.0f_d_%0.0f_x_%0.0f_z_%0.0f_subt_%0.0f_res_%0.0f_inc_%0.0fdeg_pol_%0.0fdeg%s', typestring, L(k), D(j), S(i), S(i), substrate_thickness(p), grid_size, Laser_angle(m), Laser_pol(n), endstring);
                        
                        fprintf('Working on %s\n', runningfilenamebase)
                        
                        load(runningfilenamebase)
                        u_int = u(:,1)';
                        v_int = v(1,:);
                        
                        integ_series(k) = trapz(v_int, trapz(u_int,Fh, 1));
                        length_series(k) = 1340 - 100 - L(k);
                        %length_series(k) = L(k);

                    end
                end
            end
        end
        
    end
end
disp('Done')
figure();
plot(length_series,integ_series,'k')
title('s=300 Source - Substrate Distance Background Emission Subtracted','Interpreter','latex')
%title('s=300 vs 400 - Background Emission Subtracted','Interpreter','latex')
xlabel('$$L$$ (nm)','interpreter','latex')
ylabel('Intensity (arb units)','interpreter','latex')
axis square
%legend('7 wires/um2', '11 wires/um2')