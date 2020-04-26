function [backgroundFile, numberrepeats] = backgroundFieldLogicIntegration(inputFileName, max_NWLength)
if nargin == 1   % if the number of inputs equals 1
    if strfind(inputFileName, 'maxL_2000')
        max_NWLength = 2000;
    else
        disp('Assuming a NW length of 1000 nm')
        max_NWLength = 1000;
    end
end
if max_NWLength == 2000
    if ~isempty(strfind(inputFileName, 's_300_')) || ~isempty(strfind(inputFileName, 'x_300_')) || ~isempty(strfind(inputFileName, 'x300_'))
        backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_300_maxL_2000_res_10_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_300_maxL_2000_200nmtop.mat'
        numberrepeats = 20;
     elseif~isempty(strfind(inputFileName, 's_350_')) || ~isempty(strfind(inputFileName, 'x_350_')) || ~isempty(strfind(inputFileName, 'x350_'))
        backgroundFile = 'S:\Background\3D_d40_x350_z350_maxL2000_dx0_dz0_dd0_1x1_background\background_fields_350_maxL_2000_200nmtop_1x1.mat'
        numberrepeats = 17; % This isn't quite right, so absolute intensity is slightly off, but the normalized is totally fine
    elseif~isempty(strfind(inputFileName, 's_400_')) || ~isempty(strfind(inputFileName, 'x_400_')) || ~isempty(strfind(inputFileName, 'x400_'))
        if~isempty(strfind(inputFileName, '2x2'))
                backgroundFile = 'S:\Background\3D_d40_x400_z400_maxL2000_dx0.0.0.0_dz0.0.0.0_dd0.40.0.0_2x2_background\background_fields_400_maxL_2000_200nmtop_2x2.mat'
                numberrepeats = 8;
        else
                backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_400_maxL_2000_res_10_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_400_maxL_2000_200nmtop.mat'
                numberrepeats = 15;
        end
    elseif~isempty(strfind(inputFileName, 's_500_')) || ~isempty(strfind(inputFileName, 'x_500_')) || ~isempty(strfind(inputFileName, 'x500_'))
        backgroundFile = 'S:\Background\3D_d40_x500_z500_maxL2000_dx0_dz0_dd0_1x1_background\background_fields_500_maxL_2000_200nmtop_1x1.mat'
        numberrepeats = 12;
    elseif~isempty(strfind(inputFileName, 's_600_')) || ~isempty(strfind(inputFileName, 'x_600_')) || ~isempty(strfind(inputFileName, 'x600_'))
        backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_600_maxL_2000_res_10_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_600_maxL_2000_200nmtop.mat'
        numberrepeats = 10;
    elseif~isempty(strfind(inputFileName, 's_700_')) || ~isempty(strfind(inputFileName, 'x_700_')) || ~isempty(strfind(inputFileName, 'x700_'))
        backgroundFile = 'S:\Background\3D_d40_x700_z700_maxL2000_dx0_dz0_dd0_1x1_background\background_fields_700_maxL_2000_200nmtop_1x1.mat'
        numberrepeats = 9;
    elseif~isempty(strfind(inputFileName, 's_800_')) || ~isempty(strfind(inputFileName, 'x_800_')) || ~isempty(strfind(inputFileName, 'x800_'))
        backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_x_800_z_800_maxL_2000_res_10_1x1_200nmtop_200nmsub_background\background_fields_800_maxL_2000_200nmtop.mat'
        numberrepeats = 8; % This isn't quite right, so absolute intensity is slightly off, but the normalized is totally fine
    elseif~isempty(strfind(inputFileName, 's_1000_')) || ~isempty(strfind(inputFileName, 'x_1000_')) || ~isempty(strfind(inputFileName, 'x1000_'))
        disp('Don''t have s1000 background for 2000 nm quite yet')
        return
        
        backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_1000_maxL_2000_res_10_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_1000_maxL_2000_200nmtop.mat'
        numberrepeats = 6;
    else
        disp('Didn''t find size, make sure you have a background file!')
        return
    end
    
elseif max_NWLength == 1000
    if strfind(inputFileName, '200nmtop')
        if ~isempty(strfind(inputFileName, 's_300_')) || ~isempty(strfind(inputFileName, 'x_300_'))
            backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_300_maxL_1000_res_10_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_300_200nmtop.mat'
            numberrepeats = 20;
        elseif ~isempty(strfind(inputFileName, 's_400_')) || ~isempty(strfind(inputFileName, 'x_400_'))
            if~isempty(strfind(inputFileName, '2x2'))
                backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_x_400_z_400_maxL_1000_res_10_2x2_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_400_maxL_1000_200nmtop_2x2.mat'
            else
                backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_400_maxL_1000_res_10_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_400_200nmtop.mat'
            end
            numberrepeats = 15;
        elseif ~isempty(strfind(inputFileName, 's_600_')) || ~isempty(strfind(inputFileName, 'x_600_'))
            backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_600_maxL_1000_res_10_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_600_200nmtop.mat'
            numberrepeats = 10;
        elseif ~isempty(strfind(inputFileName, 's_800_')) || ~isempty(strfind(inputFileName, 'x_800_'))
            backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_x_800_z_800_maxL_1000_res_10_1x1_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_800_maxL_1000_200nmtop.mat'
            numberrepeats = 8;% This isn't quite right, so absolute intensity is slightly off, but the normalized is totally fine
        elseif ~isempty(strfind(inputFileName, 's_1000_')) || ~isempty(strfind(inputFileName, 'x_1000_'))
            backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_1000_maxL_1000_res_10_200nmtop_200nmsub_fixedNWsize_5nm_subres4_background\background_fields_1000_200nmtop.mat'
            numberrepeats = 6;
        else
            disp('Didn''t find size, make sure you have a background file!')
            return
        end
    else
        if ~isempty(strfind(inputFileName, 's_300_')) || ~isempty(strfind(inputFileName, 'x_300_'))
            backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_300_maxL_1000_res_10_100nmtop_2nmbackground\background_fields_300.mat'
            numberrepeats = 20;
        elseif ~isempty(strfind(inputFileName, 's_400_')) || ~isempty(strfind(inputFileName, 'x_400_'))
            backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_400_maxL_1000_res_10_100nmtop_2nmbackground\background_fields_400_res2.mat'
            numberrepeats = 15;
        elseif ~isempty(strfind(inputFileName, 's_600_')) || ~isempty(strfind(inputFileName, 'x_600_'))
            backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_600_maxL_1000_res_10_100nmtop_4nmbackground\background_fields_600.mat'
            numberrepeats = 10;
        elseif ~isempty(strfind(inputFileName, 's_1000_')) || ~isempty(strfind(inputFileName, 'x_1000_'))
            backgroundFile = 'S:\Background\3D_nonuniform_xy_d_40_s_1000_maxL_1000_res_10_100nmtop_4nmbackground\background_fields_1000.mat'
            numberrepeats = 6;
        else
            disp('Didn''t find size, make sure you have a background file!')
            return
        end
    end
else
    disp('You may accidently have changed max_NWLength, only have backgrounds for 1000 and 2000 nm max lengths')
    return
end

%Only need 1 repeat for integration
if nargin == 1
    numberrepeats = 1;
end