close all;
clear all;
clc;
slice_axis = Axis.y;
slice_location = 50; %if slice_axis = Axis.y, this is nm above source (negative if you want to slice below source height
even_interpolation = true; %Always want this to be true for integration
if even_interpolation
    interpolation_gridsize = 2;
end

calculate = 'backgroundSubtracted'; %options: backgroundSubtracted, background, raw
show_figures = false;
ygrid = true; %only applies if slice_axis == Axis.y
symmetric_color = true;

k_u = pi/633; %Shouldn't need to change this, horizontal (+x) wave vector of incident light
%k_u = 0;
k_v = 0;%Shouldn't need to change this, horizontal (+z) wave vector of incident light, z into page
%Depends on slice_axis, use pi/633 for stitching along x axis (ie slicing along z or y) otherwise, 0



%Set this folder as current directory
mainfolderpath=uigetdir('','Select folder with files');
cd(mainfolderpath)
D = dir;
D = D(~ismember({D.name}, {'.', '..'}));
for m = 1:numel(D)
    if D(m).isdir
        [backgroundFile, numberrepeats] = backgroundFieldLogicIntegration(D(m).name);

        
        load(backgroundFile)
        Background_Fields = {E_background, H_background};
        source_height_background = src_array_background.intercept;
        
        
        %Even repeats for multiples of 6000 nm (ie [20 15 10 6] repeats for [300 400 600 1000])
        
        
        folderpath = fullfile(mainfolderpath,D(m).name);
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
            length_series(k) = NWlength;
            % All where they should be...
            % domainHeight(k) = obj_array(1).shape.L(2);
            % domainHeightCenter(k) = obj_array(1).shape.cb_center(2);
            % sourceHeight(k) = src_array.intercept;
            
            
            [E_compiled, H_compiled, outu, outv, domainwidth, domaindepth, k_u, k_v] = fieldSliceBackgroundSubtract(Fields, Background_Fields, k_u, k_v, source_height_background, slice_axis, slice_location, even_interpolation, interpolation_gridsize, calculate, src_array);
            
            
            E2 = E_compiled(:, :, 1).*conj(E_compiled(:, :, 1)) + E_compiled(:, :, 2).*conj(E_compiled(:, :, 2)) + E_compiled(:, :, 3).*conj(E_compiled(:, :, 3));
            H2 = H_compiled(:, :, 1).*conj(H_compiled(:, :, 1)) + H_compiled(:, :, 2).*conj(H_compiled(:, :, 2)) + H_compiled(:, :, 3).*conj(H_compiled(:, :, 3));
            
            
            integ_series(k) = fieldExtensionIntegrationAndPlot(E2, outu, outv, numberrepeats, 0, 0, domainwidth, domaindepth, slice_axis, ygrid, false, false);
            
            
            %         fprintf('Working on S\n')
            %         S = cross(E_compiled,conj(H_compiled));
            %         S_alt = cross(real(E_compiled), real(H_compiled));
            %         fprintf('Working on S_mag\n')
            %         for i=1:size(S,2)
            %             for j=1:size(S,1)
            %                 S_mag(j, i) = sqrt(S(j, i, 1)*conj(S(j, i, 1)) + S(j, i, 2)*conj(S(j, i, 2)) + S(j, i, 3)*conj(S(j, i, 3)));
            %                 S_mag_alt(j, i) = sqrt(S_alt(j, i, 1)*conj(S_alt(j, i, 1)) + S_alt(j, i, 2)*conj(S_alt(j, i, 2)) + S_alt(j, i, 3)*conj(S_alt(j, i, 3)));
            %             end
            %         end
            %
            
        end
        
        if show_figures
            fieldExtensionIntegrationAndPlot(E_compiled(:,:,3), outu, outv, numberrepeats, k_u, 0, domainwidth, domaindepth, slice_axis, ygrid, symmetric_color, show_figures);
            fieldExtensionIntegrationAndPlot(E2, outu, outv, numberrepeats, 0, 0, domainwidth, domaindepth, slice_axis, ygrid, false, show_figures);
        end
        [length_series, length_indecies] = sort(length_series);
        integ_series = integ_series(length_indecies);
        
        
        disp('Done')
        figure();
        plot(length_series, integ_series/integ_series(1), 'ko-')
        title(D(m).name,'Interpreter','latex')
        %title('s=300 vs 400 - Background Emission Subtracted','Interpreter','latex')
        xlabel('$$L$$ (nm)','interpreter','latex')
        ylabel('Intensity (arb units)','interpreter','latex')
        axis square
        %legend('7 wires/um2', '11 wires/um2')
        drawnow
        
        %Export
        save(strcat(D(m).name, 'savefile.mat'), 'length_series', 'integ_series', 'src_array', 'obj_array', 'domainwidth', 'domaindepth', 'numberrepeats')
        integFileSave = strcat(D(m).name, '_origin.txt');
        fid=fopen(integFileSave,'wt');
        fprintf(fid, [ 'Length (nm)' '\t' 'Intensity (arb units)' '\t' 'Normalized Intensity' '\r\n']);
        fprintf(fid, '%f\t%f\t%f\r\n', [length_series' integ_series' (integ_series/integ_series(1))']');
        fclose(fid);true;
        
        cd ..
    end
end

