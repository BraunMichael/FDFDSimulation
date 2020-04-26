%% SectoralCylinder
% Concrete subclass of <GenericCylinder.html |GenericCylinder|> representing a
% cylinder with a sectoral cross section.

%%% Description
% |SectoralCylinder| represents the shape of a sectoral cylinder.  Its cross
% section is a circular sector (angular portion of a disk).  The axis of the
% cylinder should be aligned with one of the axes of the Cartesian coordinate
% system.

%%% Construction
%  shape = SectoralCylinder(normal_axis, height, center, radius, theta, d_theta)
%  shape = SectoralCylinder(normal_axis, height, center, radius, theta, d_theta, dl_max)
% 
% *Input Arguments*
%
% * |normal_axis|: axis of the cylinder.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |height|: size of the cylinder along its axis.
% * |center|: center of the cylinder in the format of |[x y z]|.  For
% |normal_axis = Axis.z|, |(x, y)| is the coordinate of the center of the
% circle.
% * |radius|: radius of the circle
% * |theta|: beginning angle of the sector in radian
% * |d_theta|: angular width of the sector in radian between -2*pi and 2*pi.
% * |dl_max|: maximum grid size allowed in the cylinder.  It can be either |[dx
% dy dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned,
% |dl_max = Inf| is used.

%%% Example
%   % Create an instance of SectoralCylinder.
%   shape = SectoralCylinder(Axis.z, 100, [0 0 50], 50, pi/6, pi/3);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <SectoralCylinder.html |CircularCylinder|>, <CircularCylinder.html
% |CircularCylinder|>, <CircularShellCylinder.html
% |CircularShellCylinder|>, <EllipticCylinder.html |EllipticCylinder|>,
% <PolyognalCylinder.html |PolygonalCylinder|>, <Shape.html |Shape|>,
% <maxwell_run.html |maxwell_run|>

classdef SectoralCylinder < GenericCylinder

	properties (SetAccess = immutable)
		lsf_th  % level set function for angles
	end

	methods
        function this = SectoralCylinder(normal_axis, height, center, radius, theta, d_theta, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(height, 'real') && height > 0, '"height" should be positive.');
			chkarg(istypesizeof(center, 'real', [1, Axis.count]), ...
				'"center" should be length-%d row vector with real elements.', Axis.count);
			chkarg(istypesizeof(radius, 'real') && radius > 0, '"radius" should be positive.');
			chkarg(istypesizeof(theta, 'real'), '"theta" should be real.');
			chkarg(istypesizeof(d_theta, 'real') && (d_theta <= 2*pi || d_theta >= -2*pi), '"d_theta" should be real between -pi and pi.');
			
			thetas = sort([theta, theta + d_theta]);
			cth = mean(thetas);  % center of theta range
			cth = mod(cth + pi, 2*pi) - pi;  % -pi <= c_th < pi (range of atan2(y,x))
			sth = diff(thetas) / 2;  % semiwidth of theta range
			
			cr = radius/2;  % center of radius range
			sr = radius/2;  % semiwidth of radius range
			
			function level = lsf_th(th)
				dth = th - cth;
				dth = mod(dth + pi, 2*pi) - pi;
				level = 1 - abs(dth./sth);
			end

			[h, v, n] = cycle(normal_axis);

			% lsf2d() can handle rho = [p q] with column vectors p and q.  The
			% level set function is the one for a rectangle defined in the
			% (theta, radius) domain.
			function level = lsf2d(p, q)
				chkarg(istypeof(p, 'real'), '"p" should be array with real elements.');
				chkarg(istypeof(q, 'real'), '"q" should be array with real elements.');
				chkarg(isequal(size(p), size(q)), '"p" and "q" should have same size.');

				lc = {p - center(h), q - center(v)};  % locations in center-of-mass coordinates
				
				theta_pt = atan2(lc{Dir.v}, lc{Dir.h});  % -pi <= theta_rho < pi
				r_pt = sqrt(lc{Dir.h}.^2 + lc{Dir.v}.^2);
				zero_r_pt = (r_pt==0);  % atan2(0,0) is not well-defined; handle such cases separately
				
				dr = r_pt - cr;
				level = min(lsf_th(theta_pt), 1 - abs(dr./sr));
				level(zero_r_pt) = 0;
			end
			
			lprim = cell(1, Axis.count);
			lprim{h} = [radius * cos(thetas) + center(h), center(h)];
			lprim{v} = [radius * sin(thetas) + center(v), center(v)];
			lprim{n} = [-height height]/2 + center(n);
			
			if lsf_th(0) > 0  % sector contains +x-direction from center
				lprim{h} = [lprim{h}, center(h) + radius];
			end
			if lsf_th(pi/2) > 0  % sector contains +y-direction from center
				lprim{v} = [lprim{v}, center(v) + radius];
			end
			if lsf_th(pi) > 0  % sector contains -x-direction from center
				lprim{h} = [lprim{h}, center(h) - radius];
			end
			if lsf_th(3*pi/2) > 0  % sector contains -y-direction from center
				lprim{v} = [lprim{v}, center(v) - radius];
			end
						
			if nargin < 7  % no dl_max
				super_args = {normal_axis, @lsf2d, lprim};
			else
				super_args = {normal_axis, @lsf2d, lprim, dl_max};
			end
			
			this = this@GenericCylinder(super_args{:});
			this.lsf_th = @lsf_th;
		end
	end
end

