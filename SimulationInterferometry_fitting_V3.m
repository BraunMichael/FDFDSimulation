
%% Interferometry Deconvolution
%Michael Braun

clear; close all; format shortG; warning('off','all')
set(0,'defaultaxesfontname','arial')
set(0,'DefaultAxesFontSize',24)

%Save new text files? (0 no, 1 yes, true false also works)
Save_file = true;

%Show intermediate graphs? (0 for no, 1 for yes)
IntermediateShow = true;



%Growth Loess Smoothing factor
Growth_loess = 0;

%Derivative Loess Smoothing factor
Derivative_loess = 0;



%File Headers
header1 = 'Length (nm)';
header2 = 'Normalized Simulated Intensity';
header3 = 'Exponential Normalized Simulated Intensity';



%% Import data from text file.
addpath(genpath('E:\Michael\Stanford\Research\Data\Reflectometry'))
cd 'E:\Michael\Stanford\Research\Data\Reflectometry'
%addpath(genpath('C:\Spectre Working Folder\Reflectometry'))
%cd 'C:\Spectre Working Folder\Reflectometry'
[simFilename, folderPath] = uigetfile('*.txt;*.dat');
cd(folderPath)
[~, simFilenameOnly, ~]=fileparts(simFilename);
clc;
fprintf('Working on %s\n', simFilenameOnly);

[lengths,rawIntensities,intensities] = SimulationInterferometry_fitting_ImportFunction(simFilename);


%% Full Graph

if IntermediateShow == 1
    figure()
    plot(lengths,intensities,'k')
    ylabel('Intensity','Interpreter','latex')
    title('Normalized Interferometry Laser Deconvolution','Interpreter','latex')
    xlabel('Length (nm)','Interpreter','latex')
    axis([0 inf -inf inf])
end

Growthf=fit(lengths,intensities,'exp2');
Growthf_coeff=coeffvalues(Growthf);
a=Growthf_coeff(1);
b=Growthf_coeff(2);
c=Growthf_coeff(3);
d=Growthf_coeff(4);

Growthexpfit = @(x) a.*exp(b.*x)+c.*exp(d.*x);
Growthf_values=Growthexpfit(lengths);

Growthdiv=intensities./Growthf_values;

figure()
plot(lengths,Growthdiv,'k')
title('Growth Interferometry Deconvolved Exponential Division','Interpreter','latex')
xlabel('Length (nm)','Interpreter','latex')
ylabel('Intensity','Interpreter','latex')
axis([0 inf -inf inf])

if IntermediateShow
    figure()
    plot(lengths,intensities,'k',lengths,Growthf_values,'r')
    ylabel('Sample Signal (Volts)','Interpreter','latex')
    title('Growth Raw Interferometry Data','Interpreter','latex')
    xlabel('Length (nm)','Interpreter','latex')
    ylabel('Intensity','Interpreter','latex')
    legend('Sample','Exponential fit','Reference','Location','Best')
    axis([0 inf -inf inf])
end

Growth_smoothed = Growthdiv;
% Growth_smoothed = smooth(lengths,Growthdiv,Growth_loess,'rloess');
Deriv1_Growth_smoothed = nderiv_fornberg(1, lengths, Growth_smoothed);
Deriv2_Growth_smoothed = nderiv_fornberg(2, lengths, Growth_smoothed);

Deriv2_Growth_doublesmoothed = Deriv2_Growth_smoothed;
%Deriv2_Growth_doublesmoothed = smooth(lengths,Deriv2_Growth_smoothed,Derivative_loess,'rloess');

[zeros_1stderiv_indices, zeros_1stderiv_xvalues] = NumericalRootsFunction(lengths, Deriv1_Growth_smoothed);
[zeros_2ndderiv_indices, zeros_2ndderiv_xvalues] = NumericalRootsFunction(lengths, Deriv2_Growth_doublesmoothed);
zeros_xindices = sort([zeros_1stderiv_indices, zeros_2ndderiv_indices]);
zeros_xvalues_exact = sort([zeros_1stderiv_xvalues, zeros_2ndderiv_xvalues])';

zeros_xvalues = lengths(zeros_xindices);
zeros_yvalues = Growthdiv(zeros_xindices);

if IntermediateShow
    figure()
    plot(lengths,Growth_smoothed-1,'k',lengths,Deriv1_Growth_smoothed,'r',lengths,Deriv2_Growth_smoothed,'b')
    title('Growth Interferometry Deconvolved Exponential Division Smoothed Derivatives','Interpreter','latex')
    xlabel('Length (nm)','Interpreter','latex')
    ylabel('Intensity','Interpreter','latex')
    axis([0 inf -inf inf])
    
    figure()
    plot(lengths,Deriv2_Growth_doublesmoothed,'k',lengths,Deriv2_Growth_smoothed,'b')
    title('Growth Interferometry Deconvolved Exponential Division Smoothed Derivatives','Interpreter','latex')
    xlabel('Length (nm)','Interpreter','latex')
    ylabel('Intensity','Interpreter','latex')
    axis([0 inf -inf inf])
end

figure()
plot(lengths,Growthdiv,'k',lengths,Growth_smoothed,'r', zeros_xvalues, zeros_yvalues, 'o')
title('Growth Interferometry Deconvolved Exponential Division Smoothed','Interpreter','latex')
xlabel('Length (nm)','Interpreter','latex')
ylabel('Intensity','Interpreter','latex')
axis([0 inf -inf inf])

length_segments = (1:max(size(zeros_xvalues)));
%segmentoffset = ((length_segments(2)*zeros_xvalues(1) - zeros_xvalues(2)*length_segments(1))/(zeros_xvalues(2)-zeros_xvalues(1)));
segmentoffset = ((length_segments(2)*zeros_xvalues_exact(1) - zeros_xvalues_exact(2)*length_segments(1))/(zeros_xvalues_exact(2)-zeros_xvalues_exact(1)));
%segmentoffset2 = length_segments(1)-zeros_xvalues(1)*((length_segments(2) - length_segments(1))/(zeros_xvalues(2)-zeros_xvalues(1)));

figure()
plot(zeros_xvalues_exact, length_segments+segmentoffset, 'o')
title('Something?','Interpreter','latex')
xlabel('Length (nm)','Interpreter','latex')
ylabel('Length (arb units)','Interpreter','latex')
axis([0 inf 0 inf])

Growth_rate = nderiv_fornberg(1, zeros_xvalues_exact, length_segments+segmentoffset);
figure()
plot(zeros_xvalues_exact, Growth_rate, 'o')
title('Growth Rate vs Time','Interpreter','latex')
xlabel('Length (nm)','Interpreter','latex')
ylabel('Growth Rate (arb units)','Interpreter','latex')
axis([0 inf 0 inf])







%Save file
if Save_file==0
else
    samplefolderpath = strcat(folderPath,simFilenameOnly,'_Processed');
    mkdir(samplefolderpath)
    cd(samplefolderpath)
   
    
    Growthfilename=strcat(simFilenameOnly,'_growth.txt');
    Processedfilename=strcat(simFilenameOnly,'_processed.txt');
    
    delete(Growthfilename)
    delete(Processedfilename)
    
    

    
    
    fid=fopen(Processedfilename,'wt');
    fprintf(fid, [ 'Length (nm)' '\t' 'Exponential Divided Signal' '\t'  'Length (arb units)' '\t' '\r\n']);
    fprintf(fid, '%f\t%f\t%f\t%f\r\n', [zeros_xvalues_exact zeros_yvalues (length_segments+segmentoffset)' ]');
    fclose(fid);true;
    
    fid=fopen(Growthfilename,'wt');
    fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' 'Smoothed ExpNormDecon Signal' '\t' '1st Derivative of Smoothed' '\t' '2nd Derivative of Smoothed' '\t' 'Smoothed 2nd Derivative of Smoothed' '\r\n']);
    fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', [lengths intensities Growthdiv Growth_smoothed Deriv1_Growth_smoothed Deriv2_Growth_smoothed Deriv2_Growth_doublesmoothed]');
    fclose(fid);true;
    
    cd(folderPath)
end