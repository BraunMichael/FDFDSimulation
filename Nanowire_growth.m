clear all; close all; clear classes; clear opts; clc;

%% Simulation Setup Parameters
inspect_only = false;  % true to inspect structure before any calculation

tasks = 16; %16 %8
nodes = 2; %2 %1
total_mem = 16; %[8 16 32 48]; %Total ram, in gigabytes, scaling with spacing, just read at same index


%appendstring = '_100nmtop_4nmbackground';
%appendstring = '_100nmtop_fixedNWsize_5nm_subres4';
%appendstring = '_200nmtop_200nmsub_fixedNWsize_5nm_subres4_2NWs';
%appendstring = '_200nmtop_200nmsub_fixNW_noTip'; %_background added with createBackgroundFile
appendstring = 'pol_z_200nmtop_200nmsub'; %_background added with createBackgroundFile


timestring = {'1:00:00'}; %{'1:00:00'; '2:00:00'; '4:00:00'; '12:00:00'};
createBackgroundFile = false;
NW_exist = true;
NW_tip = true; %Simulate with gold NW catalyst?


cluster = 'Sherlock'; %Options 'Rice' or 'Sherlock'
top_source_testing = false;
UniformNWCells = false; % True makes xyz SIZE of grid cells the same in the nanowire in all dimensions, False makes xyz ASPECT RATIO of grid cells the same in the nanowires

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

superCellnum_x = 2; % number of cells in x direction
superCellnum_z = 1; % number of cells in z direction

NW_Diameter_Delta = [0 0]; %[0 40;20 0];  %must be array of size superCellnum_x (wide) by superCellnum_z ('tall'), defines deviation from D set below for wire at given position in array
NW_xpos_Delta = [0 120]; %[0 30;-120 90];  %must be array of size superCellnum_x (wide) by superCellnum_z ('tall'), defines deviation from center of subcell in x direction at given position in array
NW_zpos_Delta = [0 0]; %[0 60;-50 140]; %must be array of size superCellnum_x (wide) by superCellnum_z ('tall'), defines deviation from center of subcell in z direction at given position in array
NWspacing_x = [300 400 600]; %[300 400 600 1000]; %This is delta=0 spacing, regardless of number of wires in supercell
NWspacing_z = [300 400 600]; %These can be an array, like how S used to be, but both must match (doesn't do all possible combinations), keeps it square if an array

max_NWlength = 2000; %2000

%% Wire Parameters
substrate_exist = true;
overall_gridsize_limit = false;
%num_lengths = 99;
%L = [0 linspace(20, 1000, num_lengths)];
%L = 1000;
%L = [160 300 760 1480 1620 1760 1940]

num_lengths = 199;
L = [0 linspace(20, 2000, num_lengths)];


%num_lengths = 197;
%L = [0 linspace(40, 2000, num_lengths)];
%num_lengths = 97;
%L = [0 linspace(40, 1000, num_lengths)];


%num_lengths = 100; %20 nm steps, run again with below values to improve to 10 nm resolution if submitting in 2 batches
%L = [0 linspace(20, 2000, num_lengths)];

%num_lengths = 99; %Complement of above to mesh toghether for 10 nm resolution
%L = linspace(30, 1990, num_lengths);


%num_lengths = 51; %20 nm steps, run again with below values to improve to 10 nm resolution if submitting in 2 batches
%L = linspace(0, 1000, num_lengths);

%num_lengths = 48; %Complement of above to mesh toghether for 10 nm resolution
%L = linspace(30, 990, num_lengths);

%num_diameters = 31;
%D = linspace(40, 100, num_diameters);
D = [40 60 80];
%D = 80;


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

if max(size(NWspacing_x))>1 || max(size(NWspacing_z))>1
    if ~isequal(NWspacing_x, NWspacing_z)
        disp('Unequal x and z domain sizes are not yet supported for an array of spacings')
        return
    end
end

if ~isequal([superCellnum_z superCellnum_x], size(NW_xpos_Delta)) || ~isequal([superCellnum_z superCellnum_x], size(NW_zpos_Delta)) || ~isequal([superCellnum_z superCellnum_x], size(NW_Diameter_Delta))
    disp('The size of the superCellnum_x by superCellnum_z supercell doesn''t match the size of NW_xpos_Delta or NW_zpos_Delta or NW_Diameter_Delta arrays')
    return
end

if ((max(max(NW_xpos_Delta))+((max(D)+max(max(NW_Diameter_Delta)))/2)) > min(NWspacing_x)/2) || ((max(max(NW_zpos_Delta))+((max(D)+max(max(NW_Diameter_Delta)))/2)) > min(NWspacing_z)/2)
    disp('A nanowire in your grid might be overlapping into the cell next to it. Examine NW_xpos_Delta, NW_zpos_Delta, NW_Diameter_Delta, xNWspacing, and zNWspacing carefully!')
    return
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
startstring = sprintf('GeNW_%s',typestring);
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
if superCellnum_x > 1 || superCellnum_z > 1
    if strcmp(wire_shape, 'hexagonal')
        disp('Multiple hexagonal wires are not yet supported, must have only a single cell for hexagonal wires.')
        return
    end
end

superCellString = sprintf('_%.0fx%.0f', superCellnum_x, superCellnum_z);

if strcmp(cluster, 'Rice')
    if nodes ~= 1
        disp('Rice does terrible with more than 1 node right now')
        return
    end
    if tasks > 16
        disp('Rice does terrible with more than 1 node right now, so this also limits to 16 tasks, but it works best with 4 (right now), I promise')
        return
    end
    if max(total_mem)>96
        disp('Rice does terrible with more than 1 node right now, so limits RAM to 96? certainly less than 128 GB, but this kills it at 96 based on state normal partition amount')
        return
    end
end

folderpath=uigetdir();
addpath(genpath(folderpath))
cd(folderpath)
for bb=1:2
    if bb == 2
        NW_tip = false; %Simulate with gold NW catalyst?
        substrate_exist = false;
        overall_gridsize_limit = true; %if true, uses max_x/y/z gridsizes below on the domain, remember y is vertical (direction of wire growth)
        L = 0;
        D = 40;
        appendstring = strcat(appendstring, '_background');
    end
    for i=1:max(size(NWspacing_x))
        mem_percpu = (total_mem(i)*1024)/tasks; %in Megabytes
        
        
        if radial_growth
            if Dim_3D
                folderstring = sprintf('3D_%sx_%.0f_z_%.0f_maxL_%.0f_res_%.0f%s%s', typestring, NWspacing_x(i), NWspacing_z(i), max_NWlength, grid_size, superCellString, appendstring);
            else
                folderstring = sprintf('2D_%sx_%.0f_z_%.0f_maxL_%.0f_res_%.0f%s%s', typestring, NWspacing_x(i), NWspacing_z(i), max_NWlength, grid_size, superCellString, appendstring);
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
            switch cluster
                case 'Sherlock'
                    fprintf(fid, '#! /bin/bash\ncd /scratch/users/mbraun7/%s/RunFiles\nchmod -R +x *\ngzip -d -r /scratch/users/mbraun7/%s/\n', folderstring, folderstring);
                case 'Rice'
                    fprintf(fid, '#! /bin/bash\ncd /fasrmshare/users/mbraun7/%s/RunFiles\nchmod -R +x *\ngzip -d -r /scratch/users/mbraun7/%s/\n', folderstring, folderstring);
            end
            fclose(fid);true;
        end
        for j=1:max(size(D))
            if ~radial_growth
                midstring = sprintf('d_%0.0f_x_%0.0f_z_%0.0f_', D(j), NWspacing_x(i), NWspacing_z(i));
                if Dim_3D
                    folderstring = sprintf('3D_%s%smaxL_%.0f_res_%.0f%s%s', startstring, midstring, max_NWlength, grid_size, superCellString, appendstring);
                else
                    folderstring = sprintf('2D_%s%smaxL_%.0f_res_%.0f%s%s', startstring, midstring, max_NWlength, grid_size, superCellString, appendstring);
                end
                cd(folderpath)
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
                pythonFolderNameOnly = 'PythonFiles';
                pythonFolderString = strcat(pwd,'\',pythonFolderNameOnly);
                if 7~=exist(pythonFolderString,'dir')
                    mkdir(pythonFolderString)
                end
                addpath(genpath(pythonFolderString))
                
            end
            for m=1:max(size(Laser_angle))
                for n=1:max(size(Laser_pol))
                    endstring = sprintf('%sres_%0.0f_inc_%0.0fdeg_pol_%0.0fdeg%s%s', midstring, grid_size, Laser_angle(m), Laser_pol(n), superCellString, appendstring);
                    if radial_growth
                        masterfilename = sprintf('masterfile_%sx_%0.0f_z_%0.0f_res_%0.0f_inc_%0.0fdeg_pol_%0.0fdeg%s%s.sh', startstring, NWspacing_x(i), NWspacing_z(i), grid_size, Laser_angle(m), Laser_pol(n), superCellString, appendstring);
                    else
                        masterfilename = sprintf('masterfile_%s%s.sh', startstring, endstring);
                    end
                    fid = fopen(masterfilename,'w');
                    switch cluster
                        case 'Sherlock'
                            fprintf(fid, '#! /bin/bash\ncd /scratch/users/mbraun7/%s/RunFiles\nchmod -R +x *\ngzip -d -r /scratch/users/mbraun7/%s/\n', folderstring, folderstring);
                        case 'Rice'
                            fprintf(fid, '#! /bin/bash\ncd /farmshare/user_data/mbraun7/%s/RunFiles\nchmod -R +x *\ngzip -d -r /farmshare/user_data/mbraun7/%s/\n', folderstring, folderstring);
                    end
                    fclose(fid);true;
                    
                    pythonRestartFile(startstring, midstring, grid_size, superCellString, endstring);
                    
                    
                    
                    
                    for q=1:max(size(src_top_dist))
                        cd (fullfolderstring)
                        parfor k=1:max(size(L))
                            runningfilenamebase = sprintf('%sL_%0.0f_%s', startstring, L(k), endstring);
                            [E, H, obj_array, src_array] = nanowire_3d_func_V2(substrate_thickness, L(k), D(j), NWspacing_x(i), NWspacing_z(i), Laser_angle(m), Laser_pol(n), max_NWlength, roughness_rate, grid_size, max_wire_xz_gridsize, max_wire_y_gridsize, max_tip_gridsize, PML_cells, z_location, opts, min_tip_src_dist, src_top_dist(q), Dim_3D, solveropts, wire_shape, UniformNWCells, TFSF, inspect_only, show_solution, save_solution, runningfilenamebase, NW_tip, top_source_testing, substrate_exist, max_x_gridsize, max_y_gridsize, max_z_gridsize, overall_gridsize_limit, superCellnum_x, superCellnum_z, NW_Diameter_Delta, NW_xpos_Delta, NW_zpos_Delta, NW_exist);
                            appendInfoMatFile(runningfilenamebase, superCellnum_x, superCellnum_z, NW_Diameter_Delta, NW_xpos_Delta, NW_zpos_Delta)
                            %                         drawnow
                        end
                        cd(fullrunfilesstring)
                        for k=1:max(size(L))
                            filenamebase = sprintf('%sL_%.0f_%sres_%0.0f%s%s', startstring, L(k), midstring, grid_size, superCellString, appendstring);
                            jobname = strcat('3D_',filenamebase);
                            runningfilenamebase = sprintf('%sL_%0.0f_%s', startstring, L(k), endstring);
                            
                            pythonWatcherFileName = strcat('watcher_',filenamebase,'.py');
                            pythonPreWatcherFileName = strcat('pre',pythonWatcherFileName);
                            autopythonfilename = strcat('autorestart_',filenamebase,'.py');
                            pythonAutoRestartFile(filenamebase, autopythonfilename, pythonWatcherFileName, pythonPreWatcherFileName, pythonFolderString, pythonFolderNameOnly, folderstring, runfilesstring, cluster, timestring{i});
                            
                            fprintf('Working on %s\n', filenamebase)
                            switch cluster
                                case 'Sherlock'
                                    %Has a 30 second delay to ensure compatibility with watcher python script
                                    preamble = sprintf('#! /bin/bash\n#\n#SBATCH --job-name=%s\n#\n#SBATCH --begin=now+30\n#SBATCH --time=%s\n#SBATCH --ntasks=%.0f\n#SBATCH --nodes=%.0f\n#SBATCH --mem-per-cpu=%.0fM\n#SBATCH --mail-type=ALL\n#SBATCH --output=log%%x_%%j \ncd /scratch/users/mbraun7/%s\necho "SLURM_JOB_ID = $SLURM_JOB_ID"\ncp %s/%s/%s %s/%s/%s_q.py\npython3 %s/%s/%s $SLURM_JOB_ID\n', jobname, timestring{i}, tasks, nodes, mem_percpu, folderstring, runfilesstring, pythonFolderNameOnly, pythonWatcherFileName, runfilesstring, pythonFolderNameOnly, strcat('watcher_',filenamebase), runfilesstring, pythonFolderNameOnly, pythonPreWatcherFileName);
                                    runningfilename = sprintf('srun fd3d -i %s\n', runningfilenamebase);
                                    
                                case 'Rice' %Need to update
                                    preamble = sprintf('#! /bin/bash\n#\n#SBATCH --job-name=%s\n#\n#SBATCH --time=%s\n#SBATCH --cpus-per-task=%.0f\n#SBATCH --ntasks=1\n#SBATCH --nodes=1\n#SBATCH --mem-per-cpu=%.0fM\n#SBATCH --mail-type=ALL\n#SBATCH --output=log%%x_%%j \ncd /farmshare/user_data/mbraun7/%s\nmodule load openmpi/3.0.0\necho "SLURM_JOB_ID = $SLURM_JOB_ID"\n', jobname, timestring{i}, tasks, mem_percpu, folderstring);
                                    runningfilename = sprintf('srun --mpi=pmi2 --cpus-per-task=%.0f fd3d -i %s\n', tasks, runningfilenamebase);
                                    
                            end
                            fid=fopen(strcat(filenamebase,'.sh'),'w'); %Make it 'wt' if you want to see the line breaks in Windows, but messes with Linux
                            
                            switch cluster
                                case 'Sherlock'
                                    autopython = sprintf('cd /scratch/users/mbraun7/%s/%s\ncp log%s_$SLURM_JOB_ID temp_$SLURM_JOB_ID\npython3 %s/%s temp_$SLURM_JOB_ID\nrm temp_$SLURM_JOB_ID\n', folderstring, runfilesstring, jobname, pythonFolderNameOnly, autopythonfilename);
                                case 'Rice'
                                    autopython = sprintf('cd /farmshare/user_data/mbraun7/%s/%s\ncp log%s_$SLURM_JOB_ID temp_$SLURM_JOB_ID\npython3 %s temp_$SLURM_JOB_ID\nrm temp_$SLURM_JOB_ID\n', folderstring, runfilesstring, jobname, autopythonfilename);
                            end
                            fprintf(fid, '%s%s%s', preamble, runningfilename, autopython);
                            fclose(fid);true;
                            
                            fid = fopen(masterfilename,'a+');
                            fprintf(fid, ['sbatch ' strcat(filenamebase,'.sh') '\nsleep 1s\n']);
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
end
cd ..
disp('Done')