cd('S:\subres4\3D_nonuniform_xy_d_80_s_300_maxL_1000_res_10_100nmtop_fixedNWsize_5nm_subres4')
load('3D_nonuniform_xy_d_80_s_300_maxL_1000_res_10_100nmtop_fixedNWsize_5nm_subres4savefile.mat')


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

NWLengths = length_series;

integ_series = integ_series/integ_series(1);
minNWLengths = min(NWLengths);
maxNWLengths = max(NWLengths);
minInteg = min(integ_series);
maxInteg = max(integ_series);

max_NWLength = 1000;

Laser_angle = 60;

Movie_framerate = 10;
%Matlab_Video_Repeats = 1;
Movie_Save = true;
Saved_Video_Repeats = 1;
pauseduration = 0;


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
for j=1:max(size(NWLengths))
    ax1=subplot(1,2,1);
    cla(ax1)
    light('Position',[1 1 1],'Style','infinite')
    lighting gouraud
    hold on
    NWModelDrawing(substrateThickness, NWLengths(j), NWDiameter, NWspacing_x, NWspacing_z, Laser_angle)
    hold off
    height = max_NWLength+substrateThickness+NWRadius;
    axis(ax1, [-0.5*height 0.5*height -0.5*height 0.5*height 0 height])
    pbaspect(ax1, [1 1 1])
    view(ax1, [1 1 1]) %or view([3 3 2])
    
    ax2=subplot(1,2,2);
    if j > 1
        lengths = [lengths NWLengths(j)];
        integ = [integ integ_series(j)];
    end
    hold on
    plot(ax2,  lengths, integ,'ko-')
    plot(ax2,  lengths(end), integ(end),'ro-')
    hold off
    title(ax2, 's=300 Source - Background Emission Subtracted - 10nm NW y grid','Interpreter','latex')
    %title('s=300 vs 400 - Background Emission Subtracted','Interpreter','latex')
    xlabel(ax2, '$$L$$ (nm)','interpreter','latex')
    ylabel(ax2, 'Intensity (arb units)','interpreter','latex')
    axis(ax2, [minNWLengths maxNWLengths minInteg maxInteg])
    
    axis(ax2, 'square')
    drawnow;
    if Movie_Save
        movieframe(j)=getframe(multiframemoviefig);
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
    avifilename=strcat('Test','_SteadyStateAnimation_E2.avi');
    delete(avifilename)
    moviefile=VideoWriter(avifilename);
    moviefile.FrameRate = Movie_framerate;
    open(moviefile)
    writeVideo(moviefile,saved_movieframe);
    close(moviefile);
end


