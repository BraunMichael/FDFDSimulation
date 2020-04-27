function pythonRestartFile(startstring, endstring, typepremidstring, typepostmidstring)
pythonfilename = 'pythonrestart.py';
fid = fopen(pythonfilename,'w');
pythonFile_1 = sprintf('from subprocess import run\nimport sys\n\ndef main():\n\targs = sys.argv[1:]\n\tfor fileNumberString in args:');
pythonFile_2 = sprintf('\n\t\tif not fileNumberString.isdigit():\n\t\t\traise ValueError(''Entered value must be an int and not be the empty string'')\n\t\t\treturn\n\t\tfileNumber = int(fileNumberString)');
pythonFile_3 = sprintf('\n\t\tfileNameStart = ''%s%s''\n\t\tfileNameEnd = ''%s%s.sh''\n\t\tfileName = fileNameStart + str(fileNumber) + fileNameEnd', startstring, typepremidstring, typepostmidstring, endstring);
pythonFile_4 = sprintf('\n\t\trun("sbatch " + fileName, shell = True)\n\t\trun("sleep 1s", shell = True)\n\nif __name__ == ''__main__'':\n\tmain()');
pythonFile = strcat(pythonFile_1,pythonFile_2,pythonFile_3,pythonFile_4);
fprintf(fid, '%s', pythonFile);
fclose(fid);true;

