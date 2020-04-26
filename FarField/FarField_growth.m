clear all; close all; clear classes; clear opts; clc;

%% Simulation Setup Parameters
pathToNtoFFV8 = '/home/user/NtoFField_packageV8';
tasks = 32;
nodes = 2;
total_mem = 96; %[4 8 16 64]; %Total ram, in gigabytes, scaling with spacing, just read at same index
mem_percpu = (total_mem*1024)/tasks; %in Megabytes

%endstring = '_L_source_top_distance_90deg';
endstring = '_200nmtop_20PML';
radial_growth = false;
Dim_3D = true; %3D vs 2D?

%% Simulation Domain Parameters
grid_size = 10;  % grid cell size in nm
UniformNWCells = false; % True makes xyz SIZE of grid cells the same in the nanowire, False makes xyz NUMBER of grid cells the same in the nanowires

%% Simulation Parameters

substrate_thickness = 200; %633/2; %Wavelength/2 %Was 100 for most
%Laser_angle = 90;
Laser_angle = 60; %Angle in degrees defined as positive from substrate, 90degrees is perpendicular to substrate
Laser_pol = 0; %For TFSF: Angle defined as positive from the x-y plane, coming out of the page
%For normally, defined as positive from the +z axis towards the +x axis. +z is 0
%num_lengths = 106;

num_lengths = 101;
%num_spacings = 5;
%num_diameters = 31;
%L = linspace(50, 600, 56);
L = linspace(0, 1000, num_lengths);
%L = linspace(0, 210, num_lengths);

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
cd(folderpath)

if (tasks/nodes > 16) || ((tasks/nodes)*mem_percpu > 120000) || (nodes > 20) %> 50 on Sherlock, 20 on Rice
    disp('Too many tasks or too many tasks per node')
    return
end

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
        if ~radial_growth
            if Dim_3D
                folderstring = sprintf('3D_%sd_%.0f_s_%.0f_maxL_%.0f_res_%.0f%s', typestring, D(j), S(i), max_NWlength, grid_size, endstring);
            else
                folderstring = sprintf('2D_%sd_%.0f_s_%.0f_maxL_%.0f_res_%.0f%s', typestring, D(j), S(i), max_NWlength, grid_size, endstring);
            end
            fullfolderstring = strcat(pwd,'\',folderstring);
            if 7~=exist(fullfolderstring,'dir')
                mkdir(fullfolderstring)
            end
            addpath(genpath(fullfolderstring))
            cd (fullfolderstring)
            runfilesstring = 'RunFiles';
            fullrunfilesstring = strcat(pwd,'\',runfilesstring);
            chmodfilename = strcat(fullrunfilesstring,'\','chmodfile.sh');
            mvfilename = strcat(fullrunfilesstring,'\','mvfile.sh');
            if 7~=exist(fullrunfilesstring,'dir')
                mkdir(fullrunfilesstring)
            end
            addpath(genpath(fullrunfilesstring))
            cd(fullrunfilesstring)
            fid = fopen(chmodfilename,'w');
            fprintf(fid, '#! /bin/bash\ncd $SCRATCH/%s\nchmod -R +x *\n', folderstring);
            fclose(fid);true;
            fid = fopen(mvfilename,'w');
            fprintf(fid, '#! /bin/bash\ncd $SCRATCH/%s\n', folderstring);
            fclose(fid);true;
        end
        for m=1:max(size(Laser_angle))
            for n=1:max(size(Laser_pol))
                for p=1:max(size(substrate_thickness))
                    if radial_growth
                        masterfilename = sprintf('masterfile_GeNW_%ss_%0.0f_res_%0.0f_inc_%0.0fdeg_pol_%0.0fdeg%s.sh', typestring, S(i), grid_size, Laser_angle(m), Laser_pol(n), endstring);
                    else
                        masterfilename = sprintf('masterfile_GeNW_%sd_%0.0f_s_%0.0f_res_%0.0f_inc_%0.0fdeg_pol_%0.0fdeg%s.sh', typestring, D(j), S(i), grid_size, Laser_angle(m), Laser_pol(n), endstring);
                    end
                    fid = fopen(masterfilename,'w');
                    fprintf(fid, '#! /bin/bash\n');
                    fclose(fid);true;
                    for k=1:max(size(L))
                        filenamebase = sprintf('GeNW_%sL_%.0f_d_%0.0f_s_%0.0f_res_%0.0f%s.sh', typestring, L(k), D(j), S(i), grid_size, endstring);
                        fprintf('Working on %s\n', filenamebase)
                        jobname = sprintf('3D_%sL_%.0f_d_%0.0f_s_%0.0f_res_%0.0f%s', typestring, L(k), D(j), S(i), grid_size, endstring);
                        runningfilenamebase = sprintf('GeNW_%sL_%0.0f_d_%0.0f_x_%0.0f_z_%0.0f_subt_%0.0f_res_%0.0f_inc_%0.0fdeg_pol_%0.0fdeg%s', typestring, L(k), D(j), S(i), S(i), substrate_thickness(p), grid_size, Laser_angle(m), Laser_pol(n), endstring);
                        
                        preamble = sprintf('#! /bin/bash\n#\n#SBATCH --job-name=%s\n#\n#SBATCH --time=45:00\n#SBATCH --ntasks=%.0f\n#SBATCH --nodes=%.0f\n#SBATCH --mem-per-cpu=%.0fM\n#SBATCH --mail-type=ALL\n#SBATCH --output=log%%x \ncd $SCRATCH/%s/folder_%s\necho "SLURM_JOB_ID = $SLURM_JOB_ID"\nmodule load matlab\n', jobname, tasks, nodes, mem_percpu, folderstring, runningfilenamebase);
                        runningfilename = sprintf('matlab -r "addpath(genpath(''%s'')); N2FF_Transform; exit"', pathToNtoFFV8);
                        
                        fid=fopen(filenamebase,'w'); %Make it 'wt' if you want to see the line breaks in Windows, but messes with Linux
                        fprintf(fid, '%s%s', preamble, runningfilename);
                        fclose(fid);true;
                        
                        fid = fopen(masterfilename,'a+');
                        fprintf(fid, ['sbatch ' filenamebase '\nsleep 1s\n']);
                        fclose(fid);true;
                        
                        fid = fopen(mvfilename,'a+');
                        fprintf(fid, 'mkdir folder_%s\nmv "%s"* "folder_%s"\n', runningfilenamebase, runningfilenamebase, runningfilenamebase);
                        fclose(fid);true;
                    end
                end
            end
        end
        if ~radial_growth
            cd ..
        end
    end
    if radial_growth
        cd ..
    end
end
disp('Done')