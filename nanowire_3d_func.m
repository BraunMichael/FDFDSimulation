function [E, H, obj_array, src_array] = nanowire_3d_func_V2(substrate_t, NWlength, NWdiameter, NWspacing_x, NWspacing_z, Laser_angle, Laser_pol, max_NWlength, roughness_rate, grid_size, max_wire_xz_gridsize, max_wire_y_gridsize, max_tip_gridsize, PML_cells, z_location, opts, min_tip_src_dist, src_top_dist, Dim_3D, solveropts, wire_shape, UniformNWCells, TFSF, inspect_only, show_solution, save_solution, filenamebase, NW_tip, top_source_testing, substrate_exist, max_x_gridsize, max_y_gridsize, max_z_gridsize, overall_gridsize_limit, superCellnum_x, superCellnum_z, NW_Diameter_Delta, NW_xpos_Delta, NW_zpos_Delta, NW_exist)


%% In case of hexagonal wire
tanangle = tan(30*pi/180);
cosangle = cos(30*pi/180);
vertexarray = (NWdiameter/2).*[-tanangle -1; tanangle -1; 1/cosangle 0; tanangle 1; -tanangle 1; -1/cosangle 0];



%% Sizes
L0 = 1e-9;  % unit of length: nm
wavelength = 633;  % wavelength

Theta = Laser_angle*pi/180;
Phi = Laser_pol*pi/180;


width = NWspacing_x * superCellnum_x;
depth = NWspacing_z * superCellnum_z;
if top_source_testing
    %height = 3000;
    %height = 1790; %20PML, 200 substrate
    %source_height = 1540; %20 PML, 200 substrate
    height = 1390 + src_top_dist;
    source_height = height - src_top_dist;
else
    %height = substrate_t + max_NWlength + NWdiameter + min_tip_src_dist + src_top_dist + PML_cells*grid_size;
    height = substrate_t + max_NWlength + 40 + min_tip_src_dist + src_top_dist + PML_cells*grid_size;
    %source_height = max(height - src_top_dist, NWlength + substrate_t + NWdiameter/2 + min_tip_src_dist);
    source_height = max(height - src_top_dist, NWlength + substrate_t + 40/2 + min_tip_src_dist);
end

solveropts.filenamebase = filenamebase;
disp(solveropts.filenamebase)


%% Materials and Colors
gray = [0.5 0.5 0.5];  % [r g b]
gold = [1 215/255 0];
blue = [200/255 200/255 1];
if UniformNWCells || max_wire_xz_gridsize>=grid_size
    NW_CellAspectRatio = 1;
else
    %NW_CellAspectRatio = NWlength/NWdiameter;
    NW_CellAspectRatio = max_wire_y_gridsize/max_wire_xz_gridsize;
end


Substrate = Box([-width/2, width/2; 0, substrate_t; -depth/2, depth/2]);
Ge = Substrate;

for i = 1:superCellnum_x
    for j = 1:superCellnum_z
        cellCenter_x = (i*(width/superCellnum_x)) - 0.5*width - 0.5*NWspacing_x;
        cellCenter_z = (j*(depth/superCellnum_z)) - 0.5*depth - 0.5*NWspacing_z;
        if NW_exist && NWlength ~= 0
            if strcmp(wire_shape, 'hexagonal')
                Nanowire = PolygonalCylinder(Axis.y, NWlength, substrate_t+0.5*NWlength, vertexarray, max_wire_xz_gridsize); %Needs work if want to have adjustable supercell with hexagonal
            else
                Nanowire = CircularCylinder(Axis.y, NWlength, [cellCenter_x+NW_xpos_Delta(j,i) substrate_t+0.5*NWlength cellCenter_z+NW_zpos_Delta(j,i)], (NWdiameter+NW_Diameter_Delta(j,i))/2, max_wire_xz_gridsize.*[1 NW_CellAspectRatio 1]);
            end
            Ge = [Ge, Nanowire];
            if roughness_rate ~= 0 %This is super broken, need to add cellCenter and all kinds of things
                roughness_size = NWlength*roughness_rate;
                %Need to decide on constant number (ie just grow in size, need maxNW_length*roughness_rate below then), or if they "eat"
                %eachother (currently, using roughness_size)
                numroughperside = floor((width/2)/roughness_size);
                roughness_locations_temp = linspace(0, width/2, numroughperside+2);
                roughness_locations = roughness_locations_temp(2:(numroughperside+1));
                Roughness = [];
                for h=1:max(size(roughness_locations))
                    Roughness = horzcat(Roughness, Hemisphere([roughness_locations(h) substrate_t 0], (roughness_size/2).*[1 1 1], Axis.y));
                    Roughness = horzcat(Roughness, Hemisphere([-roughness_locations(h) substrate_t 0], (roughness_size/2).*[1 1 1], Axis.y));
                    if Dim_3D
                        for i=1:max(size(roughness_locations))
                            Roughness = horzcat(Roughness, Hemisphere([roughness_locations(h) substrate_t roughness_locations(i)], (roughness_size/2).*[1 1 1], Axis.y));
                            Roughness = horzcat(Roughness, Hemisphere([-roughness_locations(h) substrate_t roughness_locations(i)], (roughness_size/2).*[1 1 1], Axis.y));
                        end
                    end
                end
                Ge = [Ge, Roughness];
            end
        end
        if i == 1 && j == 1
            Nanowire_Tip = Hemisphere([cellCenter_x+NW_xpos_Delta(j,i) substrate_t+NWlength cellCenter_z+NW_zpos_Delta(j,i)], (NWdiameter+NW_Diameter_Delta(j,i))/2.*[1 1 1], Axis.y, max_tip_gridsize);
        else
            Nanowire_Tip = [Nanowire_Tip, Hemisphere([cellCenter_x+NW_xpos_Delta(j,i) substrate_t+NWlength cellCenter_z+NW_zpos_Delta(j,i)], (NWdiameter+NW_Diameter_Delta(j,i))/2.*[1 1 1], Axis.y, max_tip_gridsize)];
        end
    end
end

if Dim_3D
    domain_size = [-width/2, width/2; 0, height; -depth/2, depth/2];
else
    domain_size = [-width/2, width/2; 0, height; 0, grid_size];
end


if TFSF
    obj_string = 'SOBJ';
    source_string = 'SRCJ';
    source = TFSFPlaneSrc([-width/2 width/2;0 source_height;0 grid_size], [cos(Theta) -sin(Theta) 0], [sin(Theta) cos(Theta) sin(Phi)]);
else
    obj_string = 'OBJ';
    source_string = 'SRCJ';
    %source = PlaneSrc(Axis.y, source_height, Axis.x, 10, pi/2-Theta, wavelength);
    %source = PlaneSrc(Axis.x, -100, Axis.y, 10, Theta, wavelength);
    %source = PlaneSrc(Axis.y, source_height, Axis.z, 10, pi/2-Theta, wavelength);
    source = PlaneSrc(Axis.y, source_height, Phi, 10, pi/2-Theta, wavelength);
end

if overall_gridsize_limit
    domain_grid_size = [max_x_gridsize max_y_gridsize max_z_gridsize];
else
    domain_grid_size = grid_size.*[1 1 1];
end

%% Calculate Solution
if NW_tip && substrate_exist
    if strcmp(solveropts.method, 'inputfile')
        [~, ~, obj_array, src_array] = maxwell_run(...
            'OSC', L0, wavelength, ...
            'DOM', {'vacuum', 'none', 1.0}, domain_size, domain_grid_size, [BC.p BC.e BC.p], [0 PML_cells*grid_size 0],...
            obj_string, ...
            {'Palik/Ge', gray}, Ge, ...
            {'CRC/Au', gold}, Nanowire_Tip, ...
            source_string, source, solveropts, inspect_only);
        E = 0;
        H = 0;
        
    else
        [E, H, obj_array, src_array, J] = maxwell_run(...
            'OSC', L0, wavelength, ...
            'DOM', {'vacuum', 'none', 1.0}, domain_size, domain_grid_size, [BC.p BC.e BC.p], [0 PML_cells*grid_size 0],...
            obj_string, ...
            {'Palik/Ge', gray}, Ge, ...
            {'CRC/Au', gold}, Nanowire_Tip, ...
            source_string, source, solveropts, inspect_only);
    end
elseif ~NW_tip && substrate_exist
    if strcmp(solveropts.method, 'inputfile')
        [~, ~, obj_array, src_array] = maxwell_run(...
            'OSC', L0, wavelength, ...
            'DOM', {'vacuum', 'none', 1.0}, domain_size, domain_grid_size, [BC.p BC.e BC.p], [0 PML_cells*grid_size 0],...
            obj_string, ...
            {'Palik/Ge', gray}, Ge, ...
            {'vacuum', blue, 1.0}, Nanowire_Tip, ... % Didn't use to have this line, now have vacuum stand in material
            source_string, source, solveropts, inspect_only);
        E = 0;
        H = 0;
        
    else
        [E, H, obj_array, src_array, J] = maxwell_run(...
            'OSC', L0, wavelength, ...
            'DOM', {'vacuum', 'none', 1.0}, domain_size, domain_grid_size, [BC.p BC.e BC.p], [0 PML_cells*grid_size 0],...
            obj_string, ...
            {'Palik/Ge', gray}, Ge, ...
            source_string, source, solveropts, inspect_only);
    end
elseif ~NW_tip && ~substrate_exist
    if strcmp(solveropts.method, 'inputfile')
        [~, ~, obj_array, src_array] = maxwell_run(...
            'OSC', L0, wavelength, ...
            'DOM', {'vacuum', 'none', 1.0}, domain_size, domain_grid_size, [BC.p BC.e BC.p], [0 PML_cells*grid_size 0],...
            source_string, source, solveropts, inspect_only);
        E = 0;
        H = 0;
        
    else
        [E, H, obj_array, src_array, J] = maxwell_run(...
            'OSC', L0, wavelength, ...
            'DOM', {'vacuum', 'none', 1.0}, domain_size, domain_grid_size, [BC.p BC.e BC.p], [0 PML_cells*grid_size 0],...
            source_string, source, solveropts, inspect_only);
    end
end

if ~inspect_only && ~strcmp(solveropts.method, 'inputfile') && show_solution
    %% Visualize the solution.
    figure
    
    vis2d(E{Axis.x}, Axis.z, z_location, obj_array, src_array, opts)
    figure
    vis2d(E{Axis.y}, Axis.z, z_location, obj_array, src_array, opts)
    figure
    vis2d(E{Axis.z}, Axis.z, z_location, obj_array, src_array, opts)
    figure
    vis2d(H{Axis.x}, Axis.z, z_location, obj_array, src_array, opts)
    figure
    vis2d(H{Axis.y}, Axis.z, z_location, obj_array, src_array, opts)
    figure
    vis2d(H{Axis.z}, Axis.z, z_location, obj_array, src_array, opts)
    
    %Calculate energy density
    %u = 0.5 .* (PhysC.eps0.*(E*conj(E))) + (H*conj(H))./PhysC.mu0;
    %% Calculate the scatterd power
    %s_power = powerflux_box(E,H, [-sm sm; -sm sm; 0 grid_size]);
    %fprintf('scattering cross section = %e\n', s_power/grid_size);  % intensity of incident wave = 0.5
end
if save_solution && ~strcmp(solveropts.method, 'inputfile') && ~inspect_only
    filename = strcat(filenamebase, '.mat');
    save(filename, 'E', 'H', 'obj_array', 'src_array', 'J')
end
if strcmp(solveropts.method, 'inputfile') && ~inspect_only
    save(filenamebase, 'obj_array', 'src_array');
end
