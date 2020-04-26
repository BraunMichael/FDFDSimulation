clear all; close all; clear classes; clear opts; clc;

%% Simulation Setup Parameters
inspect_only = false;  % true to inspect structure before any calculation

tasks = 16;
nodes = 2;
total_mem = [8 16 32 48]; %Total ram, in gigabytes, scaling with spacing, just read at same index


%endstring = '_100nmtop_4nmbackground';
%endstring = '_100nmtop_fixedNWsize_5nm_subres4';
%endstring = '_200nmtop_200nmsub_fixedNWsize_5nm_subres4_2NWs';
endstring = '_200nmtop_200nmsub_fixedNWsize_5nm_subres4'; %_background added with createBackgroundFile

timestring = {'1:00:00'; '2:00:00'; '4:00:00'; '12:00:00'};
createBackgroundFile = true;

top_source_testing = false;
UniformNWCells = false; % True makes xyz SIZE of grid cells the same in the nanowire, False makes xyz NUMBER of grid cells the same in the nanowires

radial_growth = false;
save_solution = true;
show_solution = true;

%% Simulation Domain Parameters
grid_size = 10;  % grid cell size in nm, 10 normal
PML_cells = 10; %Number of grid cells for PML, mostly 10
min_tip_src_dist = 150; %mostly 150
%src_top_dist = 100; %mostly 50, sometimes 100
src_top_dist = 200;
max_wire_xz_gridsize = 4; %mostly 2, Utilizes non-uniform grid, minimum grid size in nm, inf uses grid_size, looks like needs to be 1 for hexagonal wire
max_wire_y_gridsize = 5; %mostly 5, 2 or 10
max_tip_gridsize = 4; %mostly 2 %Utilizes non-uniform grid, minimum grid size in nm, inf uses grid_size
num_NWs = 1; % number of NW in simulation domain
NW2_Diameter_Delta = 20;
NW2_xpos_Delta = 0;
NW2_zpos_Delta = 0;
xNWspacing = 800;
zNWspacing = 400;


S = [300 400 600 1000];
%S = [300 400];
%S = 300;

max_NWlength = 2000; %2000

if createBackgroundFile
    NW_tip = false; %Simulate with gold NW catalyst?
    substrate_exist = false;
    overall_gridsize_limit = true; %if true, uses max_x/y/z gridsizes below on the domain, remember y is vertical (direction of wire growth)
    L = 0;
    D = 40;
    endstring = strcat(endstring, '_background');
else
    NW_tip = true; %Simulate with gold NW catalyst?
    substrate_exist = true;
    overall_gridsize_limit = false;
    %num_lengths = 99;
    %L = [0 linspace(20, 1000, num_lengths)];
    %L = 1000;
    
    num_lengths = 199;
    L = [0 linspace(20, 2000, num_lengths)];
    
    %num_lengths = 97;
    %L = [0 linspace(40, 1000, num_lengths)];
    
    %num_lengths = 56;
    %L = linspace(450, 1000, num_lengths);
    
    %num_diameters = 31;
    %D = linspace(40, 100, num_diameters);
    D = [40 60 80 160];
    %D = 40;
end

%only used if createBackgroundFile is true
max_x_gridsize = 4;
max_y_gridsize = 10;
max_z_gridsize = 4;

%% Simulation Parameters

TFSF = false; %Total field scattered field source/method?
%substrate_thickness = 100; %633/2; %Wavelength/2 %Was 100 for most
substrate_thickness = 200;
Laser_angle = 60; %Angle in degrees defined as positive from substrate, 90degrees is perpendicular to substrate
Laser_pol = 0; %For TFSF: Angle defined as positive from the x-y plane, coming out of the page
%For normally, defined as positive from the +z axis towards the +x axis. +z is 0

z_location = 0;
roughness_rate = 0; %Hemisphere size as fraction of length, measured as 1/70, 0 disables roughness



solveropts.method = 'inputfile'; %'direct' or 'inputfile'
solveropts.maxit = 10e6; %Default is 1e6
solveropts.tol = 1e-6; %Default is 1e-6
solveropts.F0 = 'rand';  % 'rand' is the other choice
wire_shape = 'circular'; %Must be 'hexagonal' or 'circular'. Only works for 1 wire/cell currently
Dim_3D = true; %3D vs 2D?
opts.withobjsrc = true;
opts.withabs = false;

if max_NWlength < max(L)
    disp('Your stated max_NWlength is less than the actual maximum length (max(L))')
    return
end
if num_NWs > 1 && max(size(S)) > 1
    disp('Multiple NWs in domain doesn''t support arrays of spacing values yet!')
    return
end
if num_NWs > 1
    S = min(xNWspacing, zNWspacing);
end


if (max_wire_xz_gridsize>=grid_size) && (max_tip_gridsize>=grid_size)
    typestring = '';
elseif UniformNWCells
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
if ~xor(strcmp(solveropts.method, 'inputfile'),strcmp(solveropts.method, 'direct'))
    disp('You must choose either direct or inputfile as a solveropts.method. Check spelling and that they are strings!')
    return
end
if ~xor(strcmp(solveropts.F0, 'zero'),strcmp(solveropts.F0, 'rand'))
    disp('You must choose either zero or rand as a solveropts.F0 initialization. Check spelling and that they are strings!')
    return
end
if ~xor(strcmp(wire_shape, 'circular'),strcmp(wire_shape, 'hexagonal'))
    disp('You must choose either circular or hexagonal as the nanowire shape. Check spelling and that they are strings!')
    return
end

folderpath=uigetdir();
addpath(genpath(folderpath))
cd(folderpath)

for i=1:max(size(S))
    mem_percpu = (total_mem(i)*1024)/tasks; %in Megabytes
    
    if radial_growth
        if Dim_3D
            folderstring = sprintf('3D_%ss_%.0f_maxL_%.0f_res_%.0f%s', typestring, S(i), max_NWlength, grid_size, endstring);
        else
            folderstring = sprintf('2D_%ss_%.0f_maxL_%.0f_res_%.0f%s', typestring, S(i), max_NWlength, grid_size, endstring);
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
        if 7~=exist(fullrunfilesstring,'dir')
            mkdir(fullrunfilesstring)
        end
        addpath(genpath(fullrunfilesstring))
        cd(fullrunfilesstring)
        fid = fopen(chmodfilename,'w');
        fprintf(fid, '#! /bin/bash\ncd /scratch/users/mbraun7/%s/RunFiles\nchmod -R +x *\ngzip -d -r /scratch/users/mbraun7/%s/\n', folderstring, folderstring);
        fclose(fid);true;
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
            if 7~=exist(fullrunfilesstring,'dir')
                mkdir(fullrunfilesstring)
            end
            addpath(genpath(fullrunfilesstring))
            cd(fullrunfilesstring)
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
                    fprintf(fid, '#! /bin/bash\ncd /scratch/users/mbraun7/%s/RunFiles\nchmod +x *\ngzip -d -r /scratch/users/mbraun7/%s/\n', folderstring, folderstring);
                    fclose(fid);true;
                    
                    pythonfilename = 'pythonrestart.py';
                    fid = fopen(pythonfilename,'w');
                    pythonFile_1 = sprintf('from subprocess import run\nimport sys\n\ndef main():\n\targs = sys.argv[1:]\n\tfor fileNumberString in args:');
                    pythonFile_2 = sprintf('\n\t\tif not fileNumberString.isdigit():\n\t\t\traise ValueError(''Entered value must be an int and not be the empty string'')\n\t\t\treturn\n\t\tfileNumber = int(fileNumberString)');
                    pythonFile_3 = sprintf('\n\t\tfileNameStart = ''GeNW_%sL_''\n\t\tfileNameEnd = ''_d_%0.0f_s_%0.0f_res_%0.0f%s.sh''\n\t\tfileName = fileNameStart + str(fileNumber) + fileNameEnd',typestring, D(j), S(i), grid_size, endstring);
                    pythonFile_4 = sprintf('\n\t\trun("sbatch " + fileName, shell = True)\n\t\trun("sleep 1s", shell = True)\n\nif __name__ == ''__main__'':\n\tmain()');
                    pythonFile = strcat(pythonFile_1,pythonFile_2,pythonFile_3,pythonFile_4);
                    fprintf(fid, '%s', pythonFile);
                    fclose(fid);true;
                    
                    
                    for q=1:max(size(src_top_dist))
                        for k=1:max(size(L))
                            cd(fullrunfilesstring)
                            
                            
                            if num_NWs > 1
                                NWspacing_x = xNWspacing;
                                NWspacing_z = zNWspacing;
                            else
                                NWspacing_x = S(i);
                                NWspacing_z = S(i);
                            end
                            
                            filenamebase = sprintf('GeNW_%sL_%.0f_d_%0.0f_s_%0.0f_res_%0.0f%s.sh', typestring, L(k), D(j), NWspacing_x, grid_size, endstring);
                            jobname = sprintf('3D_%sL_%.0f_d_%0.0f_s_%0.0f_res_%0.0f%s', typestring, L(k), D(j), NWspacing_x, grid_size, endstring);
                            runningfilenamebase = sprintf('GeNW_%sL_%0.0f_d_%0.0f_x_%0.0f_z_%0.0f_subt_%0.0f_res_%0.0f_inc_%0.0fdeg_pol_%0.0fdeg%s', typestring, L(k), D(j), NWspacing_x, NWspacing_z, substrate_thickness(p), grid_size, Laser_angle(m), Laser_pol(n), endstring);
                            
                            
                            
                            
                            
                            autopythonfilename = sprintf('autorestart_GeNW_%sL_%.0f_d_%0.0f_s_%0.0f_res_%0.0f%s.py', typestring, L(k), D(j), NWspacing_x, grid_size, endstring);
                            fid = fopen(autopythonfilename,'w');
                            autopythonFile_1 = sprintf('from subprocess import run\nimport sys\n\ndef main():\n\targs = sys.argv[1:]\n\tlogFileName = args[0]');
                            autopythonFile_2 = sprintf('\n\twith open(logFileName, ''r'') as f:\n\t\tlogFile = f.read().strip()\n\t\tif logFile.find(''fd3d finished'') == -1:');
                            autopythonFile_3 = sprintf('\n\t\t\tfileName = ''GeNW_%sL_%.0f_d_%0.0f_s_%0.0f_res_%0.0f%s.sh''\n\t\t\trun("sbatch " + fileName, shell = True)',typestring, L(k), D(j), S(i), grid_size, endstring);
                            autopythonFile_4 = sprintf('\n\nif __name__ == ''__main__'':\n\tmain()');
                            autopythonFile = strcat(autopythonFile_1, autopythonFile_2, autopythonFile_3, autopythonFile_4);
                            fprintf(fid, '%s', autopythonFile);
                            fclose(fid);true;
                            
                            fprintf('Working on %s\n', filenamebase)
                            preamble = sprintf('#! /bin/bash\n#\n#SBATCH --job-name=%s\n#\n#SBATCH --time=%s\n#SBATCH --ntasks=%.0f\n#SBATCH --nodes=%.0f\n#SBATCH --mem-per-cpu=%.0fM\n#SBATCH --mail-type=ALL\n#SBATCH --output=log%%x_%%j \ncd /scratch/users/mbraun7/%s\necho "SLURM_JOB_ID = $SLURM_JOB_ID"\n', jobname, timestring{i}, tasks, nodes, mem_percpu, folderstring);
                            runningfilename = sprintf('srun fd3d -i %s\n', runningfilenamebase);
                            fid=fopen(filenamebase,'w'); %Make it 'wt' if you want to see the line breaks in Windows, but messes with Linux
                            
                            autopython = sprintf('cd /scratch/users/mbraun7/%s/%s\ncp log%s_$SLURM_JOB_ID temp_$SLURM_JOB_ID\npython3 %s temp_$SLURM_JOB_ID\nrm temp_$SLURM_JOB_ID\n', folderstring, runfilesstring, jobname, autopythonfilename);
                            fprintf(fid, '%s%s%s', preamble, runningfilename, autopython);
                            fclose(fid);true;
                            
                            fid = fopen(masterfilename,'a+');
                            fprintf(fid, ['sbatch ' filenamebase '\nsleep 1s\n']);
                            fclose(fid);true;
                            cd ..
                            
                            
                            [E, H, obj_array, src_array] = nanowire_3d_func(substrate_thickness(p), L(k), D(j), NWspacing_x, NWspacing_z, Laser_angle(m), Laser_pol(n), max_NWlength, roughness_rate, grid_size, max_wire_xz_gridsize, max_wire_y_gridsize, max_tip_gridsize, PML_cells, z_location, opts, min_tip_src_dist, src_top_dist(q), Dim_3D, solveropts, wire_shape, UniformNWCells, TFSF, inspect_only, show_solution, save_solution, runningfilenamebase, num_NWs, NW_tip, top_source_testing, substrate_exist, max_x_gridsize, max_y_gridsize, max_z_gridsize, overall_gridsize_limit, NW2_Diameter_Delta, NW2_xpos_Delta, NW2_zpos_Delta);
                            %drawnow
                            
                        end
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