function NWModelDrawing(substrateThickness, NWLength, NWDiameter, NWspacing_x, NWspacing_z, Laser_angle)

NWRadius = NWDiameter/2;
[xNW,yNW,zNW] = cylinder(NWRadius, 100);
zNW = zNW * NWLength + substrateThickness;
[xTip,yTip,zTip] = sphere(100);    
zTip(zTip<0) = 0;
r = NWRadius;
germanium = [0.5 0.5 0.5]; %[0.400000000 0.560784314 0.560784314]; %CPK coloring, otherwise [0.5 0.5 0.5]
gold = [1 215/255 0];

h = surf(r.*xTip,r.*yTip,r.*zTip + NWLength + substrateThickness);
set(h, 'edgecolor', 'none', 'facecolor', gold, 'FaceLighting','gouraud');

t = surf(xNW,yNW,zNW);
set(t, 'edgecolor', 'none', 'facecolor', germanium, 'FaceLighting', 'gouraud');

b = fill3(0.5*NWspacing_x*[-1 -1 1 1],0.5*NWspacing_z*[-1 1 1 -1],substrateThickness*[0 0 0 0],[0.5 0.5 0.5], ...  % Y faces
 0.5*NWspacing_x*[-1 -1 1 1],0.5*NWspacing_z*[-1 1 1 -1],substrateThickness*[1 1 1 1],germanium, ...  % Y faces
 0.5*NWspacing_x*[-1 -1 -1 -1],0.5*NWspacing_z*[-1 1 1 -1],substrateThickness*[0 0 1 1],germanium, ... % Z faces
 0.5*NWspacing_x*[1 1 1 1],0.5*NWspacing_z*[-1 1 1 -1],substrateThickness*[0 0 1 1],germanium, ...  % Z faces
 0.5*NWspacing_x*[-1 -1 1 1],0.5*NWspacing_z*[-1 -1 -1 -1],substrateThickness*[0 1 1 0],germanium,... % X faces
 0.5*NWspacing_x*[-1 -1 1 1],0.5*NWspacing_z*[1 1 1 1],substrateThickness*[0 1 1 0],germanium) ;  % X faces
set(b, 'FaceAlpha', 1, 'FaceLighting', 'gouraud') ;

