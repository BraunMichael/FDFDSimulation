function appendInfoMatFile(runningfilenamebase, superCellnum_x, superCellnum_z, NW_Diameter_Delta, NW_xpos_Delta, NW_zpos_Delta) 

load(strcat(runningfilenamebase, '.mat'))
save(strcat(runningfilenamebase, '.mat'))