%% Hemisphere
% Concrete subclass of <Shape.html |Shape|> representing an ellipsoid.

%%% Description
% |Hemisphere| represents the shape of a hemisphere.  The three axes of the
% hemisphere should be aligned with the axes of the Cartesian coordinate system.

%%% Construction
%  shape = Hemisphere(center, semiaxes, normalvectoraxis)
%  shape = Hemisphere(center, semiaxes, normalvectoraxis, dl_max)

% *Input Arguments*
%
% * |center|: center of the hemisphere in the format of |[x y z]|.
% * |semiaxes|: semiaxes of the hemisphere in the format of |[a b c]|. The
% flat is in the x-z plane, so the final number is the distance above the
% center it will extend
% * |normalvectoraxis|: The axis along which the normal vector of the flat
% of the sphere points. This should be |Axis.x|, |Axis.y|, or |Axis.z|
% * |dl_max|: maximum grid size allowed in the hemisphere.  It can be either |[dx
% dy dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned,
% |dl_max = Inf| is used.

%%% See Also
% <Sphere.html |Sphere|>, <CircularCylinder.html |CircularCylinder|>,
% <EllipticCylinder.html |EllipticCylinder|>, <Shape.html |Shape|>,
% <maxwell_run.html |maxwell_run|>

classdef Hemisphere < Shape
    
    methods
        function this = Hemisphere(center, semiaxes, normalvectoraxis, dl_max)
            chkarg(istypesizeof(center, 'real', [1, Axis.count]), ...
                '"center" should be length-%d row vector with real elements.', Axis.count);
            
            chkarg(istypesizeof(semiaxes, 'real', [1, Axis.count]) && all(semiaxes > 0), ...
                '"semiaxes" should be length-%d row vector with positive elements.', Axis.count);
            chkarg(istypesizeof(normalvectoraxis, 'Axis'), ...
                '"normalvectoraxis" should be instance of Axis (ie Axis.x, Axis.y, or Axis.z).');
            [h, v, n] = cycle(normalvectoraxis);
            
            bound(h,:) = [center(h) - semiaxes(h); center(h) + semiaxes(h)];
            bound(v,:) = [center(v) - semiaxes(v); center(v) + semiaxes(v)];
            bound(n,:) = [center(n); center(n) + semiaxes(n)];
            
            lprim = cell(1, Axis.count);
            for w = Axis.elems
                lprim{w} = bound(w,:);
            end
            
            % The level set function is basically 1 - norm((r-center)./semiaxes)
            % but it is vectorized, i.e., modified to handle r = [x y z] with
            % column vectors x, y, z.
            function level = lsf(x, y, z)
                chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
                chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
                chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
                chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');
                
                loc = {x, y, z};
                level = zeros(size(x));
                
                for v = Axis.elems
                    if v == normalvectoraxis
                        if loc{normalvectoraxis} < center(normalvectoraxis)
                            level = inf;
                        end
                        if loc{normalvectoraxis} == center(normalvectoraxis)
                            level = 0;
                        end
                    end
                    level = level + ((loc{v}-center(v)) ./ semiaxes(v)).^2;
                end
                level = 1 - sqrt(level);
            end
            
            if nargin < 4  % no dl_max
                super_args = {lprim, @lsf};
            else
                super_args = {lprim, @lsf, dl_max};
            end
            
            this = this@Shape(super_args{:});
        end
    end
end

