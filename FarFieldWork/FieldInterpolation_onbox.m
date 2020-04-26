function extracted_field_values = FieldInterpolation_onbox_v2(x_in, y_in, z_in)
echo on

%Different definition of coordinates between FD3D and NtoFField
xq = x_in * 1e9;
yq = z_in * 1e9;
zq = y_in * 1e9;
min_input_x = min(xq); max_input_x = max(xq);
min_input_y = min(yq); max_input_y = max(yq);
min_input_z = min(zq); max_input_z = max(zq);
fprintf('Min "z" (y) input is %.1f\n', min_input_y)
fprintf('Max "z" (y) input is %.1f\n', max_input_y)

load('saved_field_matrix')
delete('saved_field_matrix.mat')
k_u = pi/633; %Shouldn't need to change this, horizontal (+x) wave vector of incident light
k_v = 0;%Shouldn't need to change this, horizontal (+z) wave vector of incident light, z into page

Fields = {E, H};
compiled_fields = [];
for field = Fields
    for component = Axis.elems
        fprintf('Working on %s\n', field{1}{component}.name)
        primx = field{1}{component}.grid3d.lall{1,1};
        primy = field{1}{component}.grid3d.lall{2,1};
        primz = field{1}{component}.grid3d.lall{3,1};
        dualx = field{1}{component}.grid3d.lall{1,2};
        dualy = field{1}{component}.grid3d.lall{2,2};
        dualz = field{1}{component}.grid3d.lall{3,2};
        
        
        
        min_x = min(min(primx)); max_x = max(max(primx));
        %min_y = min(min(primy)); max_y = max(max(primy));
        min_z = min(min(primz)); max_z = max(max(primz));
        
        
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
        
        domainwidth = tempX(end) - tempX(1);
        domaindepth = tempZ(end) - tempZ(1);
        original_tempx = tempX;
        original_tempz = tempZ;
        
        field_values_temp = field{1}{component}.array;
        field_values_temp(1,:,:) = field_values_temp(end,:,:);
        field_values_temp(end,:,:) = field_values_temp(1,:,:).*exp(-1i*k_u*domainwidth);
        original_field_values = field_values_temp;
        
        if min_input_x < min_x
            x_blocksaway_raw = (min_input_x - min_x)/domainwidth;
            increment = -1;
            x_blocksaway_min = (sign(x_blocksaway_raw) * abs(fix(x_blocksaway_raw))) + increment; %increment cause Matlab
            temp_fields = original_field_values(1:end-1, :, :);
            for m = 1:abs(x_blocksaway_min)
                field_values_temp = cat(1,temp_fields.*exp(1i*k_u*m*domainwidth), field_values_temp);
                tempX = cat(2, (original_tempx(1:end-1) - m*domainwidth), tempX);
            end
        end
        if max_input_x > max_x
            x_blocksaway_raw = (max_input_x - max_x)/domainwidth;
            increment = 1;
            x_blocksaway_max = (sign(x_blocksaway_raw) * abs(fix(x_blocksaway_raw))) + increment; %increment cause Matlab
            temp_fields = original_field_values(2:end, :, :);
            for m = 1:x_blocksaway_max
                field_values_temp = cat(1,field_values_temp, temp_fields.*exp(-1i*k_u*m*domainwidth));
                tempX = cat(2, tempX, (original_tempx(2:end) + m*domainwidth));
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
                tempZ = cat(2, (original_tempz(1:end-1) - m*domaindepth), tempZ);
            end
        end
        if max_input_z > max_z
            z_blocksaway_raw = (max_input_z - max_z)/domaindepth;
            increment = 1;
            z_blocksaway_max = (sign(z_blocksaway_raw) * abs(fix(z_blocksaway_raw))) + increment; %increment cause Matlab
            temp_fields = intermediate_field_values(:, :, 2:end);
            for m = 1:z_blocksaway_max
                field_values_temp = cat(3, field_values_temp, temp_fields.*exp(-1i*k_v*m*domaindepth));
                tempZ = cat(2, tempZ, (original_tempz(2:end) + m*domaindepth));
            end
        end
        
        
        [X, Y, Z] = ndgrid(tempX, tempY, tempZ);
        try
            field_values = interpn(X, Y, Z, field_values_temp, xq', yq', zq', 'linear');
        catch
            save('testing')
            return
        end
        compiled_fields = [compiled_fields field_values']; %Ends up as [Ex, Ey, Ez, Hx, Hy, Hz]
    end
end
echo off
extracted_field_values = compiled_fields;


