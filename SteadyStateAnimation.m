%Steady State Animation
clear; close all; clc; format shortG; warning('off','all')
set(0,'defaultaxesfontname','arial')
set(0,'DefaultAxesFontSize',24)

%% Import data from mat file.
addpath(genpath('E:\Michael\Stanford\Research\Data\Simulation'))
cd 'E:\Michael\Stanford\Research\Data\Simulation'
%addpath(genpath('C:\Spectre Working Folder\Reflectometry'))
%cd 'C:\Spectre Working Folder\Reflectometry'
[Simfilename, folderpath] = uigetfile('*.mat');
cd(folderpath)
[~, Simfilename_only, ~]=fileparts(Simfilename);
clc;
load(Simfilename)
fprintf('Working on %s\n', Simfilename_only);

Matlab_Video_Repeats = 1;
Saved_Video_Repeats = 5;
pauseduration = 0;
Movie_Save = true;
Movie_framerate = 20;

maxamplitude = max(max(real(E2)));

% CALCULATE PHASE OF EACH FRAME
NPHI = 100;
phase = linspace(0,2*pi,NPHI);
% DRAW FRAMES
multiframemoviefig=figure();
%Maximize figure for better viewing
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
for nphi = 1 : NPHI
    % Add Phase to Field 
    f = E2.*exp(1i*phase(nphi));
    % Draw Field With Added Phase
    axhand = gca;
    fighand = pcolor(Xq, Yq, real(f));
    daspect(axhand, [1 1 1])
    set(fighand, 'EdgeColor', 'none')
    set(fighand, 'FaceColor', 'interp')
    title('$|E|^{2}$','Interpreter','latex')
    xlabel('Location (nm)','Interpreter','latex')
    ylabel('Location (nm)','Interpreter','latex')
    %colormap(axhand, b2r(4096))
    colormap(axhand, jet(4096))
    colorbar
    %caxis([-maxamplitude maxamplitude])
    caxis([0 maxamplitude])
    drawnow;
    if Movie_Save
        movieframe(nphi)=getframe(multiframemoviefig);
    end
end

for i=1:Saved_Video_Repeats
    if i>1
        saved_movieframe =  horzcat(saved_movieframe, movieframe);
    else
        saved_movieframe = movieframe;
    end
end

if Movie_Save
    %avifilename=strcat(Simfilename_only,'.avi');
    avifilename=strcat(Simfilename_only,'_SteadyStateAnimation_E2.avi');
    delete(avifilename)
    moviefile=VideoWriter(avifilename);
    moviefile.FrameRate = Movie_framerate;
    open(moviefile)
    writeVideo(moviefile,saved_movieframe);
    close(moviefile);
end