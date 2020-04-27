close all;
clear all;
clc;

set(0,'defaultaxesfontname','cmu serif')
set(0,'DefaultAxesFontSize',20)


slice_axis = Axis.z;
slice_location = 0; %if slice_axis = Axis.y, this is nm above source (negative if you want to slice below source height
even_interpolation = true; %Always want this to be true for integration
if even_interpolation
    interpolation_gridsize = 2;
end
show_figures = false;
calculate = 'backgroundSubtracted'; %options: backgroundSubtracted, background, raw
ygrid = true; %only applies if slice_axis == Axis.y
plottingfield = 'E'; %options: E2 or E

k_u = pi/633; %Shouldn't need to change this, horizontal (+x) wave vector of incident light
%k_u = 0;
k_v = 0;%Shouldn't need to change this, horizontal (+z) wave vector of incident light, z into page
%Depends on slice_axis, use pi/633 for stitching along x axis (ie slicing along z or y) otherwise, 0
[Simfilename, folderpath] = uigetfile('*.mat');
cd(folderpath)
[~, Simfilename_only, ~]=fileparts(Simfilename);
load(Simfilename)

NWLengths = length_series;
integ_series = integ_series/integ_series(1);
max_NWLength = 2000;
Laser_angle = 60;

minNWLengths = min(NWLengths);
maxNWLengths = max(NWLengths);
minInteg = min(integ_series);
maxInteg = max(integ_series);



Movie_framerate = 10;
%Matlab_Video_Repeats = 1;
Movie_Save = true;
Saved_Video_Repeats = 1;
pauseduration = 0;



overallMaxAbsValue = 0;
overallMaxValue = 0;
overallMinValue = 0;
lengths = 0;
integ = integ_series(1);
% DRAW FRAMES
multiframemoviefig=figure();
%Maximize figure for better viewing
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
light('Position',[1 1 1],'Style','infinite')
lighting gouraud


folderpath=uigetdir();
cd(folderpath)
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(folderpath, '*.E.h5');
theFiles = dir(filePattern);
[tmp, ind] = natsortfiles({theFiles.name});
theFiles = theFiles(ind);

for k = 1 : length(theFiles)
    solutionfilename = theFiles(k).name;
    [~, intermediate, ~]=fileparts(solutionfilename);
    [~, filenamebase, ~]=fileparts(intermediate);
    fprintf('Working on %s\n', filenamebase);
    load(filenamebase)
    [E, H] = read_output(filenamebase);
    Fields = {E, H};

    Compiled_fields = cell(2,1);
    Compiled_fields{1} = cell(3,1);
    Compiled_fields{2} = cell(3,1);
    
    if k == 1
        [backgroundFile, numberrepeats] = backgroundFieldLogicIntegration(filenamebase, max_NWLength);
        numberrepeats = round(numberrepeats/2);
        load(backgroundFile)
        Background_Fields = {E_background, H_background};
        source_height_background = src_array_background.intercept;
    end
    [E_compiled, H_compiled, outu, outv, domainwidth, domaindepth, k_u, k_v] = fieldSliceBackgroundSubtract(Fields, Background_Fields, k_u, k_v, source_height_background, slice_axis, slice_location, even_interpolation, interpolation_gridsize, calculate, src_array);
    switch plottingfield
        case 'E'
            field_values = E_compiled(:,:,3);
        case 'E2'
            E2 = E_compiled(:, :, 1).*conj(E_compiled(:, :, 1)) + E_compiled(:, :, 2).*conj(E_compiled(:, :, 2)) + E_compiled(:, :, 3).*conj(E_compiled(:, :, 3));
            H2 = H_compiled(:, :, 1).*conj(H_compiled(:, :, 1)) + H_compiled(:, :, 2).*conj(H_compiled(:, :, 2)) + H_compiled(:, :, 3).*conj(H_compiled(:, :, 3));
            field_values = E2;
    end
    
    field_comp = field_values;
    field_comp_temp = field_values(2:end,:);
    outu2 = outu;
    outv2 = outv;
    if numberrepeats > 1
        for m=2:numberrepeats
            switch plottingfield
                case 'E'
                    field_comp = cat(1,field_comp, field_comp_temp.*exp(-1i*k_u*(m-1)*domainwidth));
                    symmetric_color = true;
                case 'E2'
                    field_comp = cat(1,field_comp, field_comp_temp);
                    symmetric_color = false;
            end
            outu2 = cat(1, outu2, outu(2:end,:)+((m-1)*domainwidth));
            outv2 = cat(1, outv2, outv(2:end,:));
        end
    end
    if slice_axis == Axis.y && ygrid
        field_comp_temp = field_comp(:, 2:end);
        outu2_temp = outu2(:, 2:end);
        outv2_temp = outv2(:, 2:end);
        
        if numberrepeats > 1
            for n=2:numberrepeats
                field_comp = cat(2, field_comp, field_comp_temp);
                outu2 = cat(2, outu2, outu2_temp);
                outv2 = cat(2, outv2, outv2_temp + ((n-1)*domaindepth));
            end
        end
    end
    
    
    if symmetric_color
        maxabsvalue = max(max(abs(real(field_comp))));
        overallMaxAbsValue = max([maxabsvalue, overallMaxAbsValue]);
    else
        maxvalue = max(max(real(field_comp)));
        overallMaxValue = max([maxvalue, overallMaxValue]);
        minvalue = min(min(real(field_comp)));
        overallMinValue = min([minvalue, overallMinValue]);
    end
end % Find overall max values


for k = 1 : length(theFiles)
    solutionfilename = theFiles(k).name;
    [~, intermediate, ~]=fileparts(solutionfilename);
    [~, filenamebase, ~]=fileparts(intermediate);
    fprintf('Working on %s\n', filenamebase);
    load(filenamebase)
    [E, H] = read_output(filenamebase);
    Fields = {E, H};
    if strfind(filenamebase, 'L_0_')
        NWlength = 0;
    else
        for i=1:max(size(obj_array))
            if ~strcmp(obj_array(i).material.name, 'vacuum')
                if isa(obj_array(i).shape, 'CircularCylinder')
                    NWlength = obj_array(i).shape.L(2);
                    NWDiameter = obj_array(i).shape.L(1);
                    NWRadius = NWDiameter/2;
                end
            end
        end
    end
    Compiled_fields = cell(2,1);
    Compiled_fields{1} = cell(3,1);
    Compiled_fields{2} = cell(3,1);
    
    if k == 1
%         [backgroundFile, numberrepeats] = backgroundFieldLogicIntegration(filenamebase);
%         numberrepeats = round(numberrepeats/2);
%         load(backgroundFile)
%         Background_Fields = {E_background, H_background};
%         source_height_background = src_array_background.intercept;
        for i=1:max(size(obj_array))
            if ~strcmp(obj_array(i).material.name, 'vacuum')
                if isa(obj_array(i).shape, 'Box')
                    substrateThickness = obj_array(i).shape.L(2);
                    NWspacing_x = obj_array(i).shape.L(1);
                    NWspacing_z = obj_array(i).shape.L(3);
                end
                if isa(obj_array(i).shape, 'Hemisphere')
                    NWDiameter = obj_array(i).shape.L(1);
                    NWRadius = NWDiameter/2;
                end
            end
        end
    end
    [E_compiled, H_compiled, outu, outv, domainwidth, domaindepth, k_u, k_v] = fieldSliceBackgroundSubtract(Fields, Background_Fields, k_u, k_v, source_height_background, slice_axis, slice_location, even_interpolation, interpolation_gridsize, calculate, src_array);
    
    E2 = E_compiled(:, :, 1).*conj(E_compiled(:, :, 1)) + E_compiled(:, :, 2).*conj(E_compiled(:, :, 2)) + E_compiled(:, :, 3).*conj(E_compiled(:, :, 3));
    H2 = H_compiled(:, :, 1).*conj(H_compiled(:, :, 1)) + H_compiled(:, :, 2).*conj(H_compiled(:, :, 2)) + H_compiled(:, :, 3).*conj(H_compiled(:, :, 3));
    
    if show_figures
        fieldExtensionIntegrationAndPlot(E_compiled(:,:,3), outu, outv, numberrepeats, k_u, 0, domainwidth, domaindepth, slice_axis, ygrid, symmetric_color, show_figures);
        %fieldExtensionIntegrationAndPlot(E2, outu, outv, numberrepeats, 0, 0, domainwidth, domaindepth, slice_axis, ygrid, false, show_figures);
    end
    switch plottingfield
        case 'E'
            field_values = E_compiled(:,:,3);
            movieEndName = '_Growth_Ez.avi';
        case 'E2'
            field_values = E2;
            movieEndName = '_Growth_E2.avi';
    end
    
    field_comp = field_values;
    field_comp_temp = field_values(2:end,:);
    outu2 = outu;
    outv2 = outv;
    if numberrepeats > 1
        for m=2:numberrepeats
            switch plottingfield
                case 'E'
                    field_comp = cat(1,field_comp, field_comp_temp.*exp(-1i*k_u*(m-1)*domainwidth));
                    symmetric_color = true;
                case 'E2'
                    field_comp = cat(1,field_comp, field_comp_temp);
                    symmetric_color = false;
            end
            outu2 = cat(1, outu2, outu(2:end,:)+((m-1)*domainwidth));
            outv2 = cat(1, outv2, outv(2:end,:));
        end
    end
    if slice_axis == Axis.y && ygrid
        field_comp_temp = field_comp(:, 2:end);
        outu2_temp = outu2(:, 2:end);
        outv2_temp = outv2(:, 2:end);
        
        if numberrepeats > 1
            for n=2:numberrepeats
                field_comp = cat(2, field_comp, field_comp_temp);
                outu2 = cat(2, outu2, outu2_temp);
                outv2 = cat(2, outv2, outv2_temp + ((n-1)*domaindepth));
            end
        end
    end
    
    ax1=subplot(1,2,1);
    cla(ax1)
    light('Position',[1 1 1],'Style','infinite')
    lighting gouraud
    
    hold on
    % NWModelDrawing(substrateThickness, NWLengths(k), NWDiameter, NWspacing_x, NWspacing_z, Laser_angle)
    xoffset = domainwidth;
    
    for m=0:(numberrepeats-1)
        NWModelDrawing_rotated(substrateThickness, NWlength, NWDiameter, NWspacing_x, NWspacing_z, Laser_angle, m*xoffset, ax1)
    end
    
    fighand = pcolor(ax1, outu2, outv2, real(field_comp));
    set(fighand, 'EdgeColor', 'none')
    set(fighand, 'FaceColor', 'interp')
    colormap(ax1, b2r(4096))
    
    if symmetric_color
        colormap(ax1, b2r(4096))
        caxis(ax1, [-overallMaxAbsValue overallMaxAbsValue]) %For 0 = white (symmetric color axis)
    else
        colormap(ax1, hot(4096))
        caxis(ax1, [overallMinValue overallMaxValue])
    end
    hold off
    height = obj_array(1).shape.L(2);
    %axis(ax1, [-0.5*height 0.5*height 0 height -0.5*height 0.5*height ])
    axis(ax1, [-inf inf 0 height -0.5*height 0.5*height ])
    pbaspect(ax1, [1 1 1])
    axis equal
    view(ax1, [1 1 3])
    %view(ax1, [1 1 1])
    ax1.CameraUpVector = [0 1 0];
    light(ax1,'Position',[-1 1 3],'Style','infinite')
    %light(ax1,'Position',[-1 1 1.5],'Style','infinite')
    lighting gouraud
    ax1.XAxis.Visible = 'off';
    ax1.ZAxis.Visible = 'off';
    %ax1.YAxis.Visible = 'off';
    ylabel(ax1, '$$L$$ (nm)','interpreter','latex','HandleVisibility','off')
    set(gca,'NextPlot','replacechildren')
    ax2=subplot(1,2,2);
    
    if k > 1
        lengths = [lengths NWLengths(k)];
        integ = [integ integ_series(k)];
    end
    hold on
    plot(ax2, lengths, integ,'ko-')
    plot(ax2, lengths(end), integ(end),'ro-')
    hold off
    
    %title(ax2, 's=300 Source - Background Emission Subtracted - 5 nm NW y grid','Interpreter','latex')
    title(ax2, 'Simulated Reflectometry','Interpreter','latex')

    %title('s=300 vs 400 - Background Emission Subtracted','Interpreter','latex')
    xlabel(ax2, '$$L$$ (nm)','interpreter','latex')
    ylabel(ax2, 'Normalized Intensity','interpreter','latex')
    axis(ax2, [minNWLengths maxNWLengths minInteg maxInteg])
    axis(ax2, 'square')
    drawnow;
    if Movie_Save
        movieframe(k)=getframe(multiframemoviefig);
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
    avifilename=strcat(Simfilename_only, movieEndName);
    delete(avifilename)
    moviefile=VideoWriter(avifilename);
    moviefile.FrameRate = Movie_framerate;
    moviefile.Quality = 100;
    open(moviefile)
    writeVideo(moviefile,saved_movieframe);
    close(moviefile);
end
disp('Done')
