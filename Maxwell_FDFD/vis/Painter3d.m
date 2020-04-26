classdef Painter3d < handle
	% Painter3d is a class that is in charge of plotting Scalar3d based on
	% Scalar3d.
		
	properties (Access = private)
		V
		X
		Y
		Z
		maxamp  % maximum amplitude: max(abs(V(:)))
		isVpreped  % true if V is prepared
		isLpreped  % true if X, Y, Z are prepared
	end
	
	% Properties that affect V, X, Y, Z, and hence require draw2d() or draw3d()
	properties (Dependent)
		scalar3d
		grid3d
		withinterp
		withpml
		phase_angle
        withabs
		databound
		cmax
	end
	
	% Internal properties for the properties affecting V, X, Y, Z
	properties (Access = private)
		scalar3d_
		withinterp_
		withpml_
		phase_angle_
        withabs_
		cmax_
	end

	% Properties that do not affect V, X, Y, Z, but require draw().
	properties
		obj_array
		src_array
		cscale
        isopaque
        opacity
        withgrid
		withcolorbar
	end
	
	properties (SetAccess = private)
		intercept0  % default slice position
	end		
	
	methods (Access = private)
		function l_cell = lplot(this)
			% Return the locations where data are evaluated for plotting.  If
			% the data do not include the boundaries of the simulation domain (or
			% the PML interfaces for "withpml == false"), the boundary points
			% are added.
			l_cell = this.grid3d.lplot(this.scalar3d.gt_array, this.withinterp, this.withpml);
		end
		
		function l_cell = lvoxelbound(this)
			% Return the locations of boundaries of voxels drawn.  For data at
			% primary grid points, the voxel centers are in the simulation
			% domain including the boundary.  For data at dual grid points, the
			% voxel centers are in the simulation domain excluding the boundary.
			l_cell = this.grid3d.lvoxelbound(this.scalar3d.gt_array, this.withpml);
		end

		function calc_intercept0(this, V, lplot)
			absV = abs(V);
			this.intercept0 = NaN(1, Axis.count);
			for w = Axis.elems
				[u, v] = cycle(w);
				Vproj = sum(sum(absV, int(u)), int(v));
				Vproj = Vproj(:);
				[~, ind] = max(Vproj);
				this.intercept0(w) = lplot{w}(ind);
			end
		end
	end
	
	methods
		function this = Painter3d()
			this.withinterp_ = true;
			this.withpml_ = false;
			this.phase_angle_ = 0;
			this.withabs_ = false;
			this.cmax_ = NaN;

			this.obj_array = [];
			this.src_array = [];
			this.cscale = 1.0;
			this.isopaque = false;
			this.opacity = 0.7;
			this.withgrid = false;
			this.withcolorbar = false;

			this.isVpreped = false;
			this.isLpreped = false;
		end
		
		function databound = get.databound(this)
			databound = this.grid3d.bound_plot(this.withpml);
		end
				
		function scalar3d = get.scalar3d(this)
			scalar3d = this.scalar3d_;
		end
		
		function set.scalar3d(this, scalar3d)
			chkarg(istypesizeof(scalar3d, 'Scalar3d'), '"scalar3d" should be instance of Scalar3d.');
			if isempty(this.scalar3d) || ~isequal(this.scalar3d, scalar3d)  % Scalar3d is value class, so same reference means same contents
				this.isVpreped = false;
				this.isLpreped = false;
				this.scalar3d_ = scalar3d;
			end
		end
		
		function grid3d = get.grid3d(this)
			grid3d = this.scalar3d.grid3d;
		end
		
		function truth = get.withinterp(this)
			truth = this.withinterp_;
		end
		
		function set.withinterp(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			if this.withinterp ~= truth
				this.isVpreped = false;
				this.isLpreped = false;
				this.withinterp_ = truth;
			end
		end
		
		function truth = get.withpml(this)
			truth = this.withpml_;
		end
		
		function set.withpml(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			if this.withpml ~= truth
				this.isVpreped = false;
				this.isLpreped = false;
				this.withpml_ = truth;
			end
		end
		
		function phase_angle = get.phase_angle(this)
			phase_angle = this.phase_angle_;
		end
		
		function set.phase_angle(this, phase_angle)
			chkarg(istypesizeof(phase_angle, 'real'), '"phase_angle" should be real.');
			if this.phase_angle ~= phase_angle
				this.isVpreped = false;
				this.phase_angle_ = phase_angle;
			end
		end
		
		function truth = get.withabs(this)
			truth = this.withabs_;
		end
		
		function set.withabs(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			if this.withabs ~= truth
				this.isVpreped = false;
				this.withabs_ = truth;
			end
		end
		
		function cmax = get.cmax(this)
			cmax = this.cmax_;
		end

		function set.cmax(this, cmax)
			chkarg(istypesizeof(cmax, 'real') && (cmax >= 0 || isnan(cmax)), '"cmax" should be positive or NaN.');
			if isnan(cmax)
				this.isVpreped = false;
			end
			this.cmax_ = cmax;
		end
		
		function set.obj_array(this, obj_array)
			chkarg(istypesizeof(obj_array, 'EMObject', [1, 0]), '"obj_array" should be row vector of instances of EMObject.');
			this.obj_array = obj_array;
		end
		
		function set.src_array(this, src_array)
			chkarg(istypesizeof(src_array, 'Source', [1, 0]), '"src_array" should be row vector of instances of Source.');
			this.src_array = src_array;
		end
		
		function set.cscale(this, cscale)
			chkarg(istypesizeof(cscale, 'real') && cscale > 0, '"cscale" should be positive.');
			this.cscale = cscale;
		end
		
		function set.isopaque(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			this.isopaque = truth;
		end
		
		function set.opacity(this, opacity)
			chkarg(istypesizeof(opacity, 'real') && opacity >= 0 && opacity <= 1, '"opacity" should be >= 0 and <= 1.');
			this.opacity = opacity;
		end
		
		function set.withgrid(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			this.withgrid = truth;
		end
		
		function set.withcolorbar(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			this.withcolorbar = truth;
		end
		
		function prep_data(this)
			if ~this.isVpreped
				[array, lall] = this.scalar3d.data_expanded();
				[X0, Y0, Z0] = ndgrid(lall{:});

				lplot = this.lplot();
				[Xi, Yi, Zi] = ndgrid(lplot{:});
				this.V = interpn(X0, Y0, Z0, array, Xi, Yi, Zi);

				this.calc_intercept0(this.V, lplot);  % calculate before permuting V

				this.V = permute(this.V, int([Axis.y, Axis.x, Axis.z]));  % to be compatible with slice()

				this.maxamp = max(abs(this.V(:)));
				if isnan(this.cmax)
					this.cmax = this.maxamp;
				end

				if this.withabs
					this.V = abs(this.V);
				else
					this.V = real(exp(1i * this.phase_angle) .* this.V);
				end

				if ~this.withinterp
					% Pad the end boundaries.  The padded values are not drawn, so
					% they don't have to be accurate.
					this.V(end+1, :, :) = this.V(end, :, :);
					this.V(:, end+1, :) = this.V(:, end, :);
					this.V(:, :, end+1) = this.V(:, :, end);
				end
				
				this.isVpreped = true;
			end
			
			if ~this.isLpreped
				if this.withinterp
					l = this.lplot();
				else
					l = this.lvoxelbound();
				end
				[this.X, this.Y, this.Z] = meshgrid(l{:});

				this.isLpreped = true;
			end
		end
		
		function set_caxis(this, axes_handle)
			if this.withabs
				crange = this.cscale .* [0, this.cmax];
				caxis(axes_handle, crange);
				colormap(axes_handle, 'hot');
			else
				crange = this.cscale .* [-this.cmax, this.cmax];
				caxis(axes_handle, crange);
				colormap(axes_handle, b2r);
			end
		end
		
		function init_display(this, axes_handle)
			if nargin < 2  % no axes_handle
				axes_handle = gca;
			end

			this.prep_data();
			
			init_axes3d(axes_handle, this.grid3d, this.withinterp, this.withpml);
			this.set_caxis(axes_handle);
			
			% Should execute "this.set_caxis(axes_handle)" again in
			% draw_slice().

			if this.withcolorbar
				hc = colorbar;
				format_colorbar(hc, this.scalar3d);
			else
				colorbar off
			end
		end
				
		function surface_handle = draw_slice(this, axes_handle, normal_axis, intercept)
			if nargin < 2 || isempty(axes_handle)
				axes_handle = gca;
			end
			chkarg(ishandle(axes_handle), '"axes_handle" should be handle.');

			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			
			this.prep_data();

			if nargin < 4  || isempty(intercept)  % no intercept
% 				intercept = this.grid3d.center(normal_axis);
				intercept = this.intercept0(normal_axis);
			end
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');			
			chkarg(intercept >= this.databound(normal_axis, Sign.n) && intercept <= this.databound(normal_axis, Sign.p), ...
				'"intercept" should be within %s-axis bound', char(normal_axis));
									
			if normal_axis == Axis.x
				surface_handle = slice(axes_handle, this.X, this.Y, this.Z, this.V, intercept, [], []);
			elseif normal_axis == Axis.y
				surface_handle = slice(axes_handle, this.X, this.Y, this.Z, this.V, [], intercept, []);
			else
				assert(normal_axis == Axis.z);
				surface_handle = slice(axes_handle, this.X, this.Y, this.Z, this.V, [], [], intercept);
			end
			
			% Strangely, slice() reset caxis to its own default value, so caxis
			% needs to be set whenever slice() is called.
			this.set_caxis(axes_handle);
            
			if ~this.withgrid
                set(surface_handle,'EdgeColor', 'none');
			elseif this.withabs
				set(surface_handle, 'EdgeColor', 'w');
			else
				set(surface_handle, 'EdgeColor', 'k');
			end
			
			if this.isopaque
                set(surface_handle,'FaceLighting', 'phong', 'AmbientStrength', 0.5)
				if this.withinterp
					set(surface_handle, 'FaceColor', 'interp');
				end
%                light('Position',[-0.6 1 -0.2],'Style','infinite');            
			else
				if this.withinterp
					set(surface_handle, 'FaceColor', 'interp', 'FaceAlpha', 'interp');
				end
                alpha(surface_handle, this.opacity);
%                 alphamap('rampdown', this.opacity);
			end
		end
		
		function patch_handle_array = draw_objsrc(this)
			patch_handle_array = draw_objsrc(this.obj_array, this.src_array, this.grid3d, this.withinterp, this.withpml);
			
			if ~this.isopaque
				alpha(patch_handle_array, this.opacity);
			end
		end		
	end	
end

