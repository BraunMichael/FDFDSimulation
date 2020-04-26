close all;
clear all;
clc;
slice_axis = Axis.z;
slice_location = 0; %if slice_axis = Axis.y, this is nm above source (negative if you want to slice below source height
even_interpolation = true;
if even_interpolation
    interpolation_gridsize = 10;
end


save_results = false;
show_figures = true;
numberrepeats = 1;
ygrid = true; %only applies if slice_axis == Axis.y

k_x = pi/633; %Shouldn't need to change this, horizontal wave vector of incident light


%Select folder containing all the files
folderpath=uigetdir();
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(folderpath)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', folderpath);
    uiwait(warndlg(errorMessage));
    return;
end
%Set this folder as current directory
cd(folderpath);

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(folderpath, '*.E.h5');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    solutionfilename = theFiles(k).name;
    [~, intermediate, ~]=fileparts(solutionfilename);
    [~, filenamebase, ~]=fileparts(intermediate);
    fprintf('Working on %s\n', filenamebase);
    load(filenamebase)
    [E, H] = read_output(filenamebase);
    Fields = {E, H};
    %NWlength = obj_array(2).shape.L(2);
    %diameter = obj_array(2).shape.L(1);
    i=1;
    Compiled_fields = cell(2,1);
    Compiled_fields{1} = cell(3,1);
    Compiled_fields{2} = cell(3,1);
    for field = Fields
        for component = Axis.elems
            fprintf('Working on %s\n', field{1}{component}.name)
            primx = field{1}{component}.grid3d.lall{1,1};
            primy = field{1}{component}.grid3d.lall{2,1};
            primz = field{1}{component}.grid3d.lall{3,1};
            dualx = field{1}{component}.grid3d.lall{1,2};
            dualy = field{1}{component}.grid3d.lall{2,2};
            dualz = field{1}{component}.grid3d.lall{3,2};
            
            interpx_temp = sort(horzcat(primx,dualx));
            interpx = interpx_temp(interpx_temp >= min(min(primx)) & interpx_temp <= max(max(primx)));
            interpy_temp = sort(horzcat(primy,dualy));
            interpy = interpy_temp(interpy_temp >= min(min(primy)) & interpy_temp <= max(max(primy)));
            interpz_temp = sort(horzcat(primz,dualz));
            interpz = interpz_temp(interpz_temp >= min(min(primz)) & interpz_temp <= max(max(primz)));
            
            domainwidth = interpx(end)-interpx(1);
            domaindepth = interpz(end)-interpz(1);
            min_x = min(min(primx)); max_x = max(max(primx));
            min_y = min(min(primy)); max_y = max(max(primy));
            min_z = min(min(primz)); max_z = max(max(primz));
            
            
            interp = cell(1,3);
            known = cell(1,3);
            switch field{1}{component}.gt_array(Axis.x)
                case 'prim'
                    tempX = primx;
                case 'dual'
                    tempX = dualx;
            end
            switch field{1}{component}.gt_array(Axis.y)
                case 'prim'
                    tempY = primy;
                case 'dual'
                    tempY = dualy;
            end
            switch field{1}{component}.gt_array(Axis.z)
                case 'prim'
                    tempZ = primz;
                case 'dual'
                    tempZ = dualz;
            end
            
            %Use permute to switch from Maxwell_FDFD to Matlab standard of Y,X,Z array
            %indices
            %field_values_temp = permute(field{1}{component}.array, [2 1 3]);
            
            field_values_temp = field{1}{component}.array;
            k_u = k_x;
            field_values_temp(1,:,:) = field_values_temp(end,:,:);
            field_values_temp(end,:,:) = field_values_temp(1,:,:).*exp(-1i*k_u*domainwidth);
            
            
            
            %field_values_temp = field_values_temp(1:end-1,:,:);
            switch slice_axis
                case Axis.x
                    k_u = 0;
                    if slice_location < min_x || slice_location > max_x
                        disp('Your slice_location is outside the range of simulation (Axis.x)')
                        return
                    end
                    if even_interpolation
                        interp{Axis.y} = linspace(min_y, max_y, round(abs((max_y-min_y)/interpolation_gridsize))+1);
                        interp{Axis.z} = linspace(min_z, max_z, round(abs((max_z-min_z)/interpolation_gridsize))+1);
                    else
                        interp{Axis.y} = interpy;
                        interp{Axis.z} = interpz;
                    end
                    interp{Axis.x} = slice_location;
                    [inarray, slice_index] = ismembertol(slice_location, tempX);
                    
                    if inarray % interp{Axis.x} = slice_location; known{Axis.x} = slice_location;
                        [~, Yq(:, :), Zq(:, :)] = ndgrid(interp{:});
                        [Yi, Zi] = ndgrid(tempY, tempZ);
                        Compiled_fields_temp(:,:) = field_values_temp(slice_index, :, :);
                    else
                        known{Axis.y} = tempY;
                        known{Axis.z} = tempZ;
                        distance_away = sort(abs(slice_location - tempX));
                        if (distance_away(1) == distance_away(2))
                            vals = find(abs(slice_location - tempX) == distance_away(1));
                            slice_index(1) = vals(1);
                            slice_index(2) = vals(2);
                        else
                            slice_index(1) = find(abs(slice_location - tempX) == distance_away(1));
                            slice_index(2) = find(abs(slice_location - tempX) == distance_away(2));
                        end
                        slice_index = sort(slice_index);
                        known{Axis.x} = [tempX(slice_index(1)), tempX(slice_index(2))];
                        [~, Yq, Zq] = ndgrid(interp{:});
                        if length(size(Yq)) == 3
                            Yq = permute(Yq, [2 3 1]);
                            Zq = permute(Zq, [2 3 1]);
                        end
                        [Xi, Yi, Zi] = ndgrid([interp{Axis.x}, interp{Axis.x}], known{Axis.y}, known{Axis.z});
                        [X, Y, Z] = ndgrid(known{:});
                        field_values = field_values_temp(slice_index(1):slice_index(2), :, :);
                        Compiled_fields_temp = interpn(X, Y, Z, field_values, Xi, Yi, Zi, 'linear');
                        Compiled_fields_temp = permute(Compiled_fields_temp,[2 3 1]);
                        Compiled_fields_temp = Compiled_fields_temp(:,:,1);
                        if length(size(Yi)) == 3
                            Yi = permute(Yi(1,:,:), [2 3 1]);
                            Zi = permute(Zi(1,:,:), [2 3 1]);
                        end
                    end
                    Compiled_fields{i}{component} = interpn(Yi, Zi, Compiled_fields_temp, Yq, Zq, 'linear')';
                    outu = Zq'; outv = Yq';
                    
                case Axis.y
                    k_u = k_x; %Depends on slice_axis, use pi/633 for stitching along x axis (ie slicing along z or y) otherwise, 0
                    slice_location_temp = src_array.intercept + slice_location;
                    if slice_location_temp < min_y || slice_location_temp > max_y
                        disp('Your slice_location is outside the range of simulation (Axis.y)')
                        return
                    end
                    if even_interpolation
                        interp{Axis.x} = linspace(min_x, max_x, round(abs((max_x-min_x)/interpolation_gridsize))+1);
                        interp{Axis.z} = linspace(min_z, max_z, round(abs((max_z-min_z)/interpolation_gridsize))+1);
                    else
                        interp{Axis.x} = interpx;
                        interp{Axis.z} = interpz;
                    end
                    interp{Axis.y} = src_array.intercept + slice_location;
                    [inarray, slice_index] = ismembertol(slice_location_temp, tempY);
                    if inarray %interp{Axis.y} = src_array.intercept + slice_location; known{Axis.y} = slice_location_temp;
                        [Xq(:,:), ~, Zq(:, :)] = ndgrid(interp{:});
                        [Xi, Zi] = ndgrid(tempX, tempZ);
                        Compiled_fields_temp(:,:) = field_values_temp(:, slice_index, :);
                    else
                        known{Axis.x} = tempX;
                        known{Axis.z} = tempZ;
                        distance_away = sort(abs(slice_location_temp - tempY));
                        if (distance_away(1) == distance_away(2))
                            vals = find(abs(slice_location_temp - tempY) == distance_away(1));
                            slice_index(1) = vals(1);
                            slice_index(2) = vals(2);
                        else
                            slice_index(1) = find(abs(slice_location_temp - tempY) == distance_away(1));
                            slice_index(2) = find(abs(slice_location_temp - tempY) == distance_away(2));
                        end
                        slice_index = sort(slice_index);
                        known{Axis.y} = [tempY(slice_index(1)), tempY(slice_index(2))];
                        [Xq, ~, Zq] = ndgrid(interp{:});
                        if length(size(Xq)) == 3
                            Xq = permute(Xq, [1 3 2]);
                            Zq = permute(Zq, [1 3 2]);
                        end
                        [Xi, Yi, Zi] = ndgrid(known{Axis.x}, [interp{Axis.y}, interp{Axis.y}], known{Axis.z});
                        [X, Y, Z] = ndgrid(known{:});
                        field_values = field_values_temp(:, slice_index(1):slice_index(2), :);
                        Compiled_fields_temp = interpn(X, Y, Z, field_values, Xi, Yi, Zi, 'linear');
                        Compiled_fields_temp = permute(Compiled_fields_temp,[1 3 2]);
                        Compiled_fields_temp = Compiled_fields_temp(:,:,1);
                        if length(size(Xi)) == 3
                            Xi = permute(Xi(:, 1, :), [1 3 2]);
                            Zi = permute(Zi(:, 1, :), [1 3 2]);
                        end
                    end
                    Compiled_fields{i}{component} = interpn(Xi, Zi, Compiled_fields_temp, Xq, Zq, 'linear');
                    outu = Xq; outv = Zq;
                    
                case Axis.z
                    k_u = k_x; %Depends on slice_axis, use pi/633 for stitching along x axis (ie slicing along z or y) otherwise, 0
                    if slice_location < min_z || slice_location > max_z
                        disp('Your slice_location is outside the range of simulation (Axis.z)')
                        return
                    end
                    if even_interpolation
                        interp{Axis.x} = linspace(min_x, max_x, round(abs((max_x-min_x)/interpolation_gridsize))+1);
                        interp{Axis.y} = linspace(min_y, max_y, round(abs((max_y-min_y)/interpolation_gridsize))+1);
                    else
                        interp{Axis.x} = interpx;
                        interp{Axis.y} = interpy;
                    end
                    interp{Axis.z} = slice_location;
                    [inarray, slice_index] = ismembertol(slice_location, tempZ);
                    if inarray %interp{Axis.z} = slice_location; known{Axis.z} = slice_location;
                        [Xq(:, :), Yq(:, :), Zq(:, :)] = ndgrid(interp{:});
                        [Xi, Yi] = ndgrid(tempX, tempY);
                        Compiled_fields_temp(:,:) = field_values_temp(:, :, slice_index);
                    else
                        known{Axis.x} = tempX;
                        known{Axis.y} = tempY;
                        distance_away = sort(abs(slice_location - tempZ));
                        if (distance_away(1) == distance_away(2))
                            vals = find(abs(slice_location - tempZ) == distance_away(1));
                            slice_index(1) = vals(1);
                            slice_index(2) = vals(2);
                        else
                            slice_index(1) = find(abs(slice_location - tempZ) == distance_away(1));
                            slice_index(2) = find(abs(slice_location - tempZ) == distance_away(2));
                        end
                        slice_index = sort(slice_index);
                        known{Axis.z} = [tempZ(slice_index(1)), tempZ(slice_index(2))];
                        [Xq(:,:), Yq(:,:), Zq(:,:)] = ndgrid(interp{:});
                        [Xi(:,:), Yi(:,:), Zi(:,:)] = ndgrid(known{Axis.x}, known{Axis.y}, interp{Axis.z});
                        [X, Y, Z] = ndgrid(known{:});
                        field_values = field_values_temp(:, :, slice_index(1):slice_index(2));
                        Compiled_fields_temp(:,:) = interpn(X, Y, Z, field_values, Xi, Yi, Zi, 'linear');
                    end
                    
                    Compiled_fields{i}{component} = interpn(Xi, Yi, Compiled_fields_temp, Xq, Yq, 'linear');
                    outu = Xq; outv = Yq;
            end
            clear slice_location_temp vals slice_index inarray field_values field_values_temp Compiled_fields_temp X Y Z Xi Yi Zi Xq Yq Zq
        end
        i = i + 1;
    end
    
    E_compiled = cat(3, Compiled_fields{1}{Axis.x}, Compiled_fields{1}{Axis.y}, Compiled_fields{1}{Axis.z});
    H_compiled = cat(3, Compiled_fields{2}{Axis.x}, Compiled_fields{2}{Axis.y}, Compiled_fields{2}{Axis.z});
    %
    %     fprintf('Working on E2 and H2\n')
    %     for i=1:size(E_compiled,2)
    %         for j=1:size(E_compiled,1)
    %             E_mag(j, i) = real(sqrt(E_compiled(j, i, 1)*conj(E_compiled(j, i, 1)) + E_compiled(j, i, 2)*conj(E_compiled(j, i, 2)) + E_compiled(j, i, 3)*conj(E_compiled(j, i, 3))));
    %             H_mag(j, i) = real(sqrt(H_compiled(j, i, 1)*conj(H_compiled(j, i, 1)) + H_compiled(j, i, 2)*conj(H_compiled(j, i, 2)) + H_compiled(j, i, 3)*conj(H_compiled(j, i, 3))));
    %         end
    %     end
    %
    %     E2 = E_mag.^2;
    %     H2 = H_mag.^2;
    %
    %     fprintf('Working on S\n')
    %     S = cross(E_compiled,conj(H_compiled));
    %     S_alt = cross(real(E_compiled), real(H_compiled));
    %     fprintf('Working on S_mag\n')
    %     for i=1:size(S,2)
    %         for j=1:size(S,1)
    %             S_mag(j, i) = sqrt(S(j, i, 1)*conj(S(j, i, 1)) + S(j, i, 2)*conj(S(j, i, 2)) + S(j, i, 3)*conj(S(j, i, 3)));
    %             S_mag_alt(j, i) = sqrt(S_alt(j, i, 1)*conj(S_alt(j, i, 1)) + S_alt(j, i, 2)*conj(S_alt(j, i, 2)) + S_alt(j, i, 3)*conj(S_alt(j, i, 3)));
    %         end
    %     end
    %
    %     if save_results
    %         fprintf('Working on saving\n')
    %         savename = 'test.mat'; %filenamebase
    %         save(savename, 'E', 'H', 'obj_array', 'src_array', 'E_compiled', 'H_compiled', 'E2', 'H2', 'S', 'S_mag')
    %     end
end

fprintf('Done, working on figures now if flag is set\n')
if show_figures
    E_comp = E_compiled(:,:,3);
    E_comp_temp = E_compiled(2:end,:,3);
    outu2 = outu;
    outv2 = outv;
    if numberrepeats > 1
        for m=2:numberrepeats
            E_comp = cat(1,E_comp, E_comp_temp.*exp(-1i*k_u*(m-1)*domainwidth));
            outu2 = cat(1, outu2, outu(2:end,:)+((m-1)*domainwidth));
            outv2 = cat(1, outv2, outv(2:end,:));
        end
    end
    if slice_axis == Axis.y && ygrid
        E_comp_temp = E_comp(:, 2:end);
        outu2_temp = outu2(:, 2:end);
        outv2_temp = outv2(:, 2:end);
        
        if numberrepeats > 1
            for n=2:numberrepeats
                E_comp = cat(2, E_comp, E_comp_temp);
                outu2 = cat(2, outu2, outu2_temp);
                outv2 = cat(2, outv2, outv2_temp + ((n-1)*domaindepth));
            end
        end
    end
    maxvalue = max(max(abs(real(E_comp))));
    
    figure
    hold on
    axhand = gca;
    %fighand = pcolor(outx, outy, E2);
    fighand = pcolor(outu2, outv2, real(E_comp));
    daspect(axhand, [1 1 1])
    set(fighand, 'EdgeColor', 'none')
    set(fighand, 'FaceColor', 'interp')
    xlabel('Location (nm)','Interpreter','latex')
    ylabel('Location (nm)','Interpreter','latex')
    colormap(axhand, b2r(4096))
    colorbar
    caxis([-maxvalue maxvalue])
    hold off
%     
%     figure
%     axis([min(min(outu2)) max(max(outu2)) min(min(outv2)) max(max(outv2))])
%     line([min(min(outu2)),max(max(outu2))], [src_array.intercept, src_array.intercept], 'Color', 'green')
%     line([min(min(outu2)),max(max(outu2))], [obj_array(3).shape.L(2), obj_array(3).shape.L(2)], 'Color', 'black')
%     for n=1:numberrepeats
%         try
%             pos = [(obj_array(2).shape.bound(1,1)+ ((n-1)*domainwidth)) obj_array(2).shape.bound(2,1) obj_array(2).shape.bound(1,2)-obj_array(2).shape.bound(1,1) obj_array(2).shape.bound(2,2)-obj_array(2).shape.bound(2,1)]; %four-element vector of the form [x y w h]
%             rectangle('Position', pos)
%         catch
%         end
%         try
%             thet = linspace(0,pi);
%             x = (obj_array(4).shape.cb_center(1) + ((n-1)*domainwidth)) + 0.5*obj_array(4).shape.L(1)*cos(thet);
%             y = obj_array(4).shape.bound(2,1) + obj_array(4).shape.L(2)*sin(thet);
%             plot(x,y,'Color',[1 215/255 0])
%         catch
%         end
%         hold off
%         drawnow
%     end
%     figpos = get(gca, 'Position'); %[xmin ymin width height]
%     figratio = figpos(3)/figpos(4);
%     set(gcf,'PaperSize',[10*figratio 10]); %[width height]
end
    %Export
    %save('fieldvis.mat', 'outu2', 'outv2', 'E_comp', 'maxvalue', 'src_array', 'obj_array', 'domainwidth', 'numberrepeats')