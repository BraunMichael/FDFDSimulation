function pythonAutoRestartFile(filenamebase, autopythonfilename, pythonWatcherFileName, pythonPreWatcherFileName, pythonFolderString, pythonFolderNameOnly, folderstring, runfilesstring, cluster, timestring)

[~, pythonWatcherScriptNameOnly, ~] = fileparts(pythonWatcherFileName);


% Need sprintf for the tabs/new lines
fid = fopen(strcat(pythonFolderString,'/',autopythonfilename),'w');
autopythonFile_1 = sprintf('from subprocess import run\nimport sys\n\ndef main():\n\targs = sys.argv[1:]\n\tlogFileName = args[0]');
autopythonFile_2 = sprintf('\n\twith open(logFileName, ''r'') as f:\n\t\tlogFile = f.read().strip()\n\t\tif logFile.find(''fd3d finished'') == -1:');
autopythonFile_3 = sprintf('\n\t\t\tfileName = ''%s.sh''\n\t\t\trun("sbatch " + fileName, shell = True)', filenamebase);
autopythonFile_4 = sprintf('\n\nif __name__ == ''__main__'':\n\tmain()');
autopythonFile = strcat(autopythonFile_1, autopythonFile_2, autopythonFile_3, autopythonFile_4);
fprintf(fid, '%s', autopythonFile);
fclose(fid);true;


fid = fopen(strcat(pythonFolderString,'/',pythonPreWatcherFileName), 'w');
pythonPreWatcherFile_1 = sprintf('from subprocess import run\nimport sys\nimport fileinput\n\ndef main():\n\targs = sys.argv[1:]\n\tif len(args) > 1:');
pythonPreWatcherFile_2 = sprintf('\n\t\traise ValueError(''Will only operate on 1 job'')\n\t\treturn\n\tjobNumber = args[0]\n\tif not jobNumber.isdigit():\n\t\traise ValueError(''Entered value must be an int and not be the empty string'')\n\t\treturn');
pythonPreWatcherFile_3 = sprintf('\n\tfor line in fileinput.input("%s/%s/%s_q.py",inplace=True):\n\t\tprint(line.replace("REPLACE_ME",jobNumber),end="")\n\tfileName = ''%s/%s/%s.sh''\n\trun("sbatch " + fileName, shell = True)\n\texit()', runfilesstring, pythonFolderNameOnly, pythonWatcherScriptNameOnly, runfilesstring, pythonFolderNameOnly, pythonWatcherScriptNameOnly);
pythonPreWatcherFile_4 = sprintf('\n\nif __name__ == ''__main__'':\n\tmain()');
pythonPreWatcherFile = strcat(pythonPreWatcherFile_1,pythonPreWatcherFile_2,pythonPreWatcherFile_3,pythonPreWatcherFile_4);
fprintf(fid, '%s', pythonPreWatcherFile);
fclose(fid);true;

fid = fopen(strcat(pythonFolderString,'/',pythonWatcherFileName), 'w');
pythonWatcherFile_1 = sprintf('import sys\nimport subprocess\n\ndef get_sec(time_str):\n\tsplitTimeString = time_str.split('':'')\n\tif len(splitTimeString) != 3:\n\t\treturn float("inf")');
pythonWatcherFile_2 = sprintf('\n\tfor timeUnit in splitTimeString:\n\t\tif not timeUnit.isdigit():\n\t\t\treturn float("inf")\n\th, m, s = time_str.split('':'')\n\treturn int(h) * 3600 + int(m) * 60 + int(s)');
pythonWatcherFile_3 = sprintf('\n\ndef main():\n\tjobNumber = "REPLACE_ME"\n\tfileName = ''%s.sh''\n\tif not jobNumber.isdigit():\n\t\traise ValueError(''Entered value must be an int and not be the empty string'')\n\t\treturn\n\thistoryOutput = subprocess.check_output("sacct --format=state,elapsed --jobs=" + jobNumber,shell=True).decode("utf-8")', filenamebase);
pythonWatcherFile_4 = sprintf('\n\tif historyOutput.find(''TIMEOUT'') != -1:\n\t\tsubprocess.run("sbatch " + fileName, shell = True)\n\t\treturn\n\tif historyOutput.find(''CANCELLED'') != -1:\n\t\thistoryOutputSplit = historyOutput.split()');
pythonWatcherFile_5 = sprintf('\n\t\tfor historyColumn in historyOutputSplit:\n\t\t\tif get_sec(historyColumn) < 120:\n\t\t\t\tsubprocess.run("sbatch " + fileName, shell = True)\n\t\t\t\treturn');
pythonWatcherFile_6 = sprintf('\n\nif __name__ == ''__main__'':\n\tmain()');
pythonWatcherFile = strcat(pythonWatcherFile_1,pythonWatcherFile_2,pythonWatcherFile_3,pythonWatcherFile_4,pythonWatcherFile_5,pythonWatcherFile_6);
fprintf(fid, '%s', pythonWatcherFile);
fclose(fid);true;


hms = strsplit(timestring,':');

jobWaitSeconds = 3600 * str2double(hms{1}) + 60 * str2double(hms{2}) + str2double(hms{3});
delaySeconds = 120;
fid = fopen(strcat(pythonFolderString,'/',pythonWatcherScriptNameOnly,'.sh'), 'w');
switch cluster
    case 'Sherlock'
        preamble = sprintf('#! /bin/bash\n#\n#SBATCH --job-name=%s\n#\n#SBATCH --begin=now+%.0f\n#SBATCH --time=00:05:00\n#SBATCH --ntasks=1\n#SBATCH --nodes=1\n#SBATCH --mem-per-cpu=2048M\n#SBATCH --mail-type=ALL\n#SBATCH --output=log%%x_%%j \ncd /scratch/users/mbraun7/%s/%s\necho "SLURM_JOB_ID = $SLURM_JOB_ID"\npython3 %s/%s_q.py\n', pythonWatcherScriptNameOnly, jobWaitSeconds+delaySeconds, folderstring, runfilesstring, pythonFolderNameOnly, pythonWatcherScriptNameOnly);
        
    case 'Rice' %Need to update
        preamble = sprintf('#! /bin/bash\n#\n#SBATCH --job-name=%s\n#\n#SBATCH --begin=now+%.0f\n#SBATCH --time=00:05:00\n#SBATCH --ntasks=1\n#SBATCH --nodes=1\n#SBATCH --mem-per-cpu=2048M\n#SBATCH --mail-type=ALL\n#SBATCH --output=log%%x_%%j \ncd /farmshare/user_data/mbraun7/%s/%s\nmodule load openmpi/3.0.0\necho "SLURM_JOB_ID = $SLURM_JOB_ID"\npython3 %s/%s_q.py\n', pythonWatcherScriptNameOnly, jobWaitSeconds+delaySeconds, folderstring, runfilesstring, pythonFolderNameOnly, pythonWatcherScriptNameOnly);
end
fprintf(fid, '%s', preamble);
fclose(fid);true;