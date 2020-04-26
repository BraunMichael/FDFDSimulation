clear all; close all; clear classes; clear opts; clc;

tasks = 32;
nodes = 2;
total_mem = 64; %Total ram, in gigabytes
    mem_percpu = (total_mem*1024)/tasks; %in Megabytes

folderpath=uigetdir();
addpath(genpath(folderpath))
cd(folderpath)

typefilename = sprintf('3D_N2FFtransform');
folderstring = sprintf('%s', typefilename);
remotefolderpath_cd = sprintf('cd $HOME/nodetesting/%s/', folderstring);
mpath_cd = sprintf('cd $HOME/NtoFField_packageV8/');


chmodfilename = strcat(folderpath,'\','chmodfile.sh');
fid = fopen(chmodfilename,'w');
fprintf(fid, '#! /bin/bash\n');
fclose(fid);true;

fullfolderstring = strcat(pwd,'\',folderstring);
if 7~=exist(fullfolderstring,'dir')
    mkdir(fullfolderstring)
end
addpath(genpath(fullfolderstring))
cd (fullfolderstring)

masterfilename = sprintf('nodetest_masterfile_GeNW_%s.sh', typefilename);
fid = fopen(masterfilename,'w');
fprintf(fid, '#! /bin/bash\n');
fclose(fid);true;

fid = fopen(chmodfilename,'a+');
fprintf(fid, '%s\nchmod +x *\n./%s', remotefolderpath_cd, masterfilename);
fclose(fid);true;

for task = tasks
    %nodes = 2.^(0:1:log2(task));
    for node = nodes
        if (task/node > 16) || ((task/node)*mem_percpu > 120000) || (node > 20) %> 50 on Sherlock, 20 on Rice
            continue
        else
            filenamebase = sprintf('GeNW_%s_%.0ftasks_%.0fnodes.sh', typefilename, task, node);
            fprintf('Working on %s\n', filenamebase)
            jobname = sprintf('%s_%.0ftasks_%.0fnodes_360', typefilename, task, node);
            preamble = sprintf('#! /bin/bash\n#\n#SBATCH --job-name=%s\n#\n#SBATCH --time=4:00:00\n#SBATCH --ntasks=%.0f\n#SBATCH --nodes=%.0f\n#SBATCH --mem-per-cpu=%.0fM\n#SBATCH --mail-type=ALL\n#SBATCH --output=log%%x \n%s\necho "SLURM_JOB_ID = $SLURM_JOB_ID"\nmodule load matlab\n', jobname, task, node, mem_percpu, mpath_cd);

            runningfilename = sprintf('matlab -r ''N2FF_Transform; exit''\n');
            
            fid=fopen(filenamebase,'w'); %Make it 'wt' if you want to see the line breaks in Windows, but messes with Linux
            fprintf(fid, '%s%s', preamble, runningfilename);
            fclose(fid);true;
            
            fid = fopen(masterfilename,'a+');
            fprintf(fid, ['sbatch ' filenamebase '\nsleep 1s\n']);
            fclose(fid);true;
            
        end
    end
end
cd ..