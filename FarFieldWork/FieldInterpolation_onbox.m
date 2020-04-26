function extracted_field_values = FieldInterpolation_onbox_v4(x_in, y_in, z_in)
echo on

%% Convert coordinate definition between FD3D and NtoFField
xq = x_in * 1e9;
yq = z_in * 1e9;
zq = y_in * 1e9;
min_input_x = min(xq); max_input_x = max(xq);
min_input_y = min(yq); max_input_y = max(yq);
min_input_z = min(zq); max_input_z = max(zq);
fprintf('Min "z" (y) input (bottom cut) is %.1f\n', min_input_y)
fprintf('Max "z" (y) input (top cut) is %.1f\n', max_input_y)

k_u = pi/633; %Shouldn't need to change this, horizontal (+x) wave vector of incident light
k_v = 0;%Shouldn't need to change this, horizontal (+z) wave vector of incident light, z into page

% Fetch background intensity
%load('background_fields_300.mat')
%load('background_fields_400.mat')
load('background_fields_300_res4.mat');
Background_Fields = {E_background, H_background};
% Fetch main intensity
load('saved_field_matrix')
delete('saved_field_matrix.mat')
Fields = {E, H};
compiled_fields = [];

for i = 1:2
    for component = Axis.elems
        if i == 1
            field = Fields{1};
            background_field = Background_Fields{1};
        elseif i == 2
            field = Fields{2};
            background_field = Background_Fields{2};
        else
            disp('Fix your loop indexing in FieldInterpolation_onbox')
            return
        end
        fprintf('Working on %s\n', field{component}.name)
        primx = field{component}.grid3d.lall{1,1};
        primy = field{component}.grid3d.lall{2,1};
        primz = field{component}.grid3d.lall{3,1};
        dualx = field{component}.grid3d.lall{1,2};
        dualy = field{component}.grid3d.lall{2,2};
        dualz = field{component}.grid3d.lall{3,2};
        
        primx_background = background_field{component}.grid3d.lall{1,1};
        primy_background = background_field{component}.grid3d.lall{2,1};
        primz_background = background_field{component}.grid3d.lall{3,1};
        dualx_background = background_field{component}.grid3d.lall{1,2};
        dualy_background = background_field{component}.grid3d.lall{2,2};
        dualz_background = background_field{component}.grid3d.lall{3,2};
        
        switch field{component}.gt_array(Axis.x)
            case 'prim'
                tempX = primx;
            case 'dual'
                tempX = dualx;
        end
        switch field{component}.gt_array(Axis.y)
            case 'prim'
                tempY = primy;
            case 'dual'
                tempY = dualy;
        end
        switch field{component}.gt_array(Axis.z)
            case 'prim'
                tempZ = primz;
            case 'dual'
                tempZ = dualz;
        end
        
        switch background_field{component}.gt_array(Axis.x)
            case 'prim'
                tempX_background = primx_background;
            case 'dual'
                tempX_background = dualx_background;
        end
        switch background_field{component}.gt_array(Axis.y)
            case 'prim'
                tempY_background = primy_background;
            case 'dual'
                tempY_background = dualy_background;
        end
        switch background_field{component}.gt_array(Axis.z)
            case 'prim'
                tempZ_background = primz_background;
            case 'dual'
                tempZ_background = dualz_background;
        end
        
        min_x = min(min(primx)); max_x = max(max(primx));
        min_y = min(min(primy)); max_y = max(max(primy));
        min_z = min(min(primz)); max_z = max(max(primz));
        
        
        
        domainwidth_background_initial = tempX_background(end) - tempX_background(1);
        domainwidth_initial = tempX(end) - tempX(1);
        domaindepth_initial = tempZ(end) - tempZ(1);
        
        if ((domainwidth_initial < 0.9*domainwidth_background_initial) || (domainwidth_initial > 1.1*domainwidth_background_initial))
            disp('Your background and main domain widths are probably different sizes, be sure to include matching domain files! If you are sure you have the right background file, comment out the return below this line!')
            return
        end
        
        field_values_temp_presubtraction = field{component}.array;
        field_values_temp_presubtraction(1,:,:) = field_values_temp_presubtraction(end,:,:);
        field_values_temp_presubtraction(end,:,:) = field_values_temp_presubtraction(1,:,:).*exp(-1i*k_u*domainwidth_initial);
        
        background_field_values_temp = background_field{component}.array;
        background_field_values_temp(1,:,:) = background_field_values_temp(end,:,:);
        background_field_values_temp(end,:,:) = background_field_values_temp(1,:,:).*exp(-1i*k_u*domainwidth_background_initial);

        

        
        tempX_clipped = tempX;
        tempX_clipped(tempX_clipped <= min_x) = [];
        tempX_clipped(tempX_clipped >= max_x) = [];
        if tempX_clipped(1) < tempX_clipped(end)
            tempX_clipped = [min_x, tempX_clipped, max_x];
        else
            tempX_clipped = [max_x, tempX_clipped, min_x];
        end
        
        
       tempY_clipped = tempY;
       tempY_clipped(tempY_clipped <= min_y) = [];
       tempY_clipped(tempY_clipped >= max_y) = [];
       if tempY_clipped(1) < tempY_clipped(end)
           tempY_clipped = [min_y, tempY_clipped, max_y];
       else
           tempY_clipped = [max_y, tempY_clipped, min_y];
       end
        
        tempZ_clipped = tempZ;
        tempZ_clipped(tempZ_clipped <= min_z) = [];
        tempZ_clipped(tempZ_clipped >= max_z) = [];
        if tempZ_clipped(1) < tempZ_clipped(end)
            tempZ_clipped = [min_x, tempZ_clipped, max_x];
        else
            tempZ_clipped = [max_x, tempZ_clipped, min_x];
        end
        
                
        domainwidth = abs(max_x - min_x);
        domaindepth = abs(max_z - min_z);        
        

        
        [X_background, Y_background, Z_background] = ndgrid(tempX_background, tempY_background, tempZ_background);
        [X, Y, Z] = ndgrid(tempX, tempY, tempZ);
        [X_clipped, Y_clipped, Z_clipped] = ndgrid(tempX_clipped, tempY_clipped, tempZ_clipped);
        
        
        try
            background_field_values = interpn(X_background, Y_background, Z_background, background_field_values_temp, X_clipped, Y_clipped, Z_clipped, 'linear');
        catch
            disp('Broke on background_field_values')
            i
            component
            save('testing_background_field_values.mat', '-v7.3')
            return
        end
        
        try
            field_values_temp_presubtraction_clipped = interpn(X, Y, Z, field_values_temp_presubtraction, X_clipped, Y_clipped, Z_clipped, 'linear');
        catch
            disp('Broke on field_values_temp clipping')
            i
            component
            save('testing_background_field_values.mat', '-v7.3')
            return
        end
        
        %% Subtract background
        background_field_values(Y_clipped < source_height_background) = 0;
        field_values_temp = field_values_temp_presubtraction_clipped - background_field_values;
        %save('bg_fieldsubtest.mat', 'field_values_temp', 'field_values_temp_presubtraction','background_field_values','source_height_background','domainwidth','domainwidth_background_initial', 'tempX_background', 'tempY_background', 'tempZ_background', 'tempX', 'tempY', 'tempZ', 'min_x', 'min_z', 'tempX_clipped', 'tempY_clipped', 'tempZ_clipped')
        
        %% Extend domain
        original_tempx = tempX_clipped;
        original_tempz = tempZ_clipped;
        original_field_values = field_values_temp;
        if min_input_x < min_x
            x_blocksaway_raw = (min_input_x - min_x)/domainwidth;
            increment = -1;
            x_blocksaway_min = (sign(x_blocksaway_raw) * abs(fix(x_blocksaway_raw))) + increment; %increment cause Matlab
            temp_fields = original_field_values(1:end-1, :, :);
            for m = 1:abs(x_blocksaway_min)
                field_values_temp = cat(1,temp_fields.*exp(1i*k_u*m*domainwidth), field_values_temp);
                tempX_clipped = cat(2, (original_tempx(1:end-1) - m*domainwidth), tempX_clipped);
            end
        end
%        save('testing_field_values.mat', '-v7.3')
%        return
        if max_input_x > max_x
            x_blocksaway_raw = (max_input_x - max_x)/domainwidth;
            increment = 1;
            x_blocksaway_max = (sign(x_blocksaway_raw) * abs(fix(x_blocksaway_raw))) + increment; %increment cause Matlab
            temp_fields = original_field_values(2:end, :, :);
            for m = 1:x_blocksaway_max
                field_values_temp = cat(1,field_values_temp, temp_fields.*exp(-1i*k_u*m*domainwidth));
                tempX_clipped = cat(2, tempX_clipped, (original_tempx(2:end) + m*domainwidth));
            end
        end
        intermediate_field_values = field_values_temp;
        if min_input_z < min_z
            z_blocksaway_raw = (min_input_z - min_z)/domaindepth;
            increment = -1;
            z_blocksaway_min = (sign(z_blocksaway_raw) * abs(fix(z_blocksaway_raw))) + increment; %increment cause Matlab
            temp_fields = intermediate_field_values(:, :, 1:end-1);
            for m = 1:abs(z_blocksaway_min)
                field_values_temp = cat(3,temp_fields.*exp(1i*k_v*m*domaindepth), field_values_temp);
                tempZ_clipped = cat(2, (original_tempz(1:end-1) - m*domaindepth), tempZ_clipped);
            end
        end
        if max_input_z > max_z
            z_blocksaway_raw = (max_input_z - max_z)/domaindepth;
            increment = 1;
            z_blocksaway_max = (sign(z_blocksaway_raw) * abs(fix(z_blocksaway_raw))) + increment; %increment cause Matlab
            temp_fields = intermediate_field_values(:, :, 2:end);
            for m = 1:z_blocksaway_max
                field_values_temp = cat(3, field_values_temp, temp_fields.*exp(-1i*k_v*m*domaindepth));
                tempZ_clipped = cat(2, tempZ_clipped, (original_tempz(2:end) + m*domaindepth));
            end
        end
        
        [X, Y, Z] = ndgrid(tempX_clipped, tempY_clipped, tempZ_clipped);
        try
            field_values = interpn(X, Y, Z, field_values_temp, xq', yq', zq', 'linear');
        catch
            disp('Broke on field_values')
            i
            component
            save('testing_field_values.mat', '-v7.3')
            return
        end
        compiled_fields = [compiled_fields field_values']; %Ends up as [Ex, Ey, Ez, Hx, Hy, Hz]
    end
end


echo off
extracted_field_values = compiled_fields;


