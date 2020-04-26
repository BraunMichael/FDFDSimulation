clear all; close all; clear classes; clc;

%% Set flags.
isnew = false;
inspect_only = true;
srcfilenamebase = 'coax';

%% Solve the system.
if isnew
	%%
	gray = [0.5 0.5 0.5];  % [r g b]
	solveropts.method = 'inputfile';
	solveropts.filenamebase = mfilename;

	srcfilename = srcfilenamebase;
	load(srcfilename, 'src_array');
	modalsrc = src_array;

	dL = 20;
	dl = 2;
	xn = -600 - 10*dL;
	xp = 300 + 10*dL;
	yn = xn;
	yp = xp;
	zn = -300 - 10*dL;
	zp = 300 + 10*dL;
	wo = 150;
	wi = 50;
	t = 25;
	T = 50;
	wvlen = 1550;
	[E, H, obj_array, src_array, J] = maxwell_run(...
		'OSC', 1e-9, wvlen, ...
		'DOM', {'vacuum', 'none', 1.0}, [xn xp; yn yp; zn zp], dL, BC.m, 10*dL, ...
		'OBJ', ...
			{'vacuum', 'none', 1.0}, ...
				Box([xn xp; yn yp; -wo/2-T wo/2+t], [dL dL dL/2]), ...
				Plane(Axis.z, 0), ...
				Plane(Axis.x, -400), ...
				Plane(Axis.y, -400), ...
			{['CRC', filesep, 'Ag'], gray}, ...
				Box([xn xp; yn yp; zn -wo/2]), ...
				Box([-wo/2-t wo/2+t; yn wo/2+t; -wo/2 wo/2+t], [dl dL dl]), ...
				Box([xn wo/2+t; -wo/2-t wo/2+t; -wo/2 wo/2+t], [dL dl dl]), ...
			{['Palik', filesep, 'Si'], 'm'}, ...
				Box([xn xp; yn yp; zn -wo/2-T]), ...
			{['Palik', filesep, 'SiO2'], 'y'}, ...
				Box([-wo/2 wo/2; yn wo/2; -wo/2 wo/2], [dl dL dl]), ...
				Box([xn wo/2; -wo/2 wo/2; -wo/2 wo/2], [dL dl dl]), ...
			{['CRC', filesep, 'Ag'], gray}, ...
				Box([-wi/2 wi/2; yn wi/2; -wi/2 wi/2], [dl dL dl]), ...
				Box([xn wi/2; -wi/2 wi/2; -wi/2 wi/2], [dL dl dl]), ...
		'SRCJ', modalsrc, ...
		solveropts, inspect_only);

	if ~inspect_only
		save(solveropts.filenamebase, 'J', 'obj_array', 'src_array');
	end
else
	load(mfilename);
	[E, H] = read_output(mfilename);
end

%%
if ~inspect_only
	figure;
	clear opts;
	opts.withabs = false;
	opts.isopaque = false;
% 	visall(H{Axis.z}, obj_array, src_array, opts);
	vis3d(H{Axis.z}, [], 0, 0, opts);
	axis(gca, [-200 200 -200 200 -200 200]);
	colormap(gca, 'jet');
end
