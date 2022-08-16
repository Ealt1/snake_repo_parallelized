classdef Rec<handle

	properties (SetAccess=protected)
		dim;	% Dimension of rec, scalar
		xmin;	% Row vector of lower bounds
		xmax;	% Row vector of upper bounds
		ap;		% List of atomic propositions in cell
	end

	methods
		function r = Rec(X, ap)
			% Rec: Create a hyperrectangle
			% 
			% SYNTAX
			% ------
			%
			%	rec = Rec(X);
			% 	rec = Rec(X,ap);
			%
			% INPUT
			% -----
			%	
			%	X 	n x 2 matrix defining min and max along each dimension.
			%		Class: double
			%	ap 	(array of) APs of of the hyperrectangle. Default: []
			%		Class: int
			%
			% PROPERTIES
			% ------
			%
			% 	dim 	dimension
			%	xmin 	array (1 x dim) of lower bounds for each dimension
			%	xmax 	array (1 x dim) of upper bounds for each dimension
			% 
			% SEE ALSO
			% -------
			% 
			% add_ap, remove_ap, isFullDim, isEmptySet, getFlatDims, getFullDims, getMidpoint, 
			% volume, split, isInside, contains, contains_strictly, getVertices, getVertexI,
			% getFacetI, isNeighbor, mldivide, plot, projection, toPoly
			if nargin<2
				ap = [];
			end

			if size(X,1) ~= 2
				X = X';
			end
			if size(X,1) ~= 2
				error('Wrong dimensions of X')
			end
			r.dim = size(X,2);
			r.xmin = X(1,:);
			r.xmax = X(2,:);
			r.ap = ap;
		end

		function add_ap(rec, ap)
			% Add an atomic proposition
			rec.ap = union(rec.ap, ap);
		end

		function remove_ap(rec, ap)
			% Remove an atomic proposition
			rec.ap = setdiff(rec.ap, ap);
		end

		function ret = isFullDim(rec)
			% Returns true if Rec is fully dimensional
			ret = all(rec.xmin<rec.xmax);
		end

		function ind = getFlatDims(rec)
			% Returns indices of "flat" dimensionssuch that x_{i,min} = x_{i,max}
			if rec.isEmptySet
				ind = [];
				return;
			end
			ind = find(rec.xmin == rec.xmax);
		end

		function ind = getFullDims(rec)
			% Returns indices of "full" dimensions i such that x_{i,min} < x_{i,max}
			if rec.isEmptySet
				ind = [];
				return;
			end
			ind = find(rec.xmin < rec.xmax);
		end

		function ret = isEmptySet(rec)
			% Returns true if Rec` is the empty set
            % if rec is an array of Rec's, return 1 if the union is empty 
            ret = [];
            for i = 1:1:length(rec)
                ret = [ret, any(rec(i).xmin>rec(i).xmax)];
            end
            ret = all(ret);
		end

		function mid = getMidpoint(recs)
			% Returns the middle point of a Rec 
            mid = []; 
            for i = 1:1:length(recs)
                rec = recs(i); 
                mid = [mid; (rec.xmin + rec.xmax)/2];
            end
            mid = mid'; 
		end

		function vol = volume(rec)
			% Computes volume
			if length(rec)>1
				vol = zeros(1,length(rec));
				for i = 1:length(rec)
					vol(i) = volume(rec(i));
				end
				return;
			end
			if isEmptySet(rec)
                vol = 0;
            else
                vol = prod(rec.xmax-rec.xmin);
            end
		end

		function [part1, part2] = split(rec, dim)
			% Split the Rec along dimension dim
			xmin = rec.xmin;
			xmax = rec.xmax;
			midval = (xmax(dim) + xmin(dim))/2;
			xminspl = xmin; xminspl(dim) = midval;
			xmaxspl = xmax; xmaxspl(dim) = midval;
			part1 = Rec([xmin; xmaxspl], rec.ap);
			part2 = Rec([xminspl; xmax], rec.ap);
		end

		function ret = eq(rec1, rec2)
			% Overloaded == operator
			if length(rec1)>1 || length(rec2)>1
				ret1 = mldivide(rec1, rec2);
				ret2 = mldivide(rec2, rec1);
				ret = isempty(ret1) && isempty(ret2);
				return;
			end
			ret = all(rec1.xmin == rec2.xmin) && all(rec1.xmax == rec2.xmax);
		end

		function ret = ne(rec1, rec2)
			% Overloaded ~= operator
			ret = ~eq(rec1,rec2);
		end

		function ret = isInside(rec, point)
			% Return true if point is inside rec
			if length(point) ~= rec.dim
				error('isInside: dimension mismatch')
			end
			point = reshape(point,1,length(point));
			ret = all(rec.xmin<=point) && all(point<=rec.xmax);
		end

		function ret = contains(rec1, rec2)
			% returns true if rec2 is contained in rec1
			ret = all(rec2.xmax<=rec1.xmax) && all(rec2.xmin>=rec1.xmin);
        end
        
        function ret = multiple_contains_one(recs, rec)
            % assume all full dimensional
            epsilon = volume(rec)/100^rec.dim;
            Vol_intersection = 0;
            for i = 1:1:length(recs)
                Vol_intersection = Vol_intersection ...
                    + volume(intersect(recs(i), rec));
            end
            ret = ((volume(rec) - Vol_intersection)<=epsilon); 
        end

		function ret = contains_strictly(rec1, rec2)
		    % Return true if rec2 does not touch the boundaries of rec1 along dimension dim
		    ret = all(rec1.xmax < rec2.xmax) && all(rec1.xmin>rec2.xmin);
        end
        
        function collapse(rec, cdim)
            % inputs: cdim, row vector
            % assume: cdim is subset of 1:1:rec.dim
            ncdim = setdiff(1:1:rec.dim, cdim);
            rec.xmin = rec.xmin(ncdim);
            rec.xmax = rec.xmax(ncdim);
            rec.dim = rec.dim - length(cdim);
        end

        function vert = getVertices(rec)
            % Returns a matrix where each line is the coordinate of a vertex of rec

            if rec.isEmptySet
                vert = [];
                return;
            end

            if rec.dim == 1
                vert = [rec.xmin rec.xmax];
                return;
            end

            if rec.dim == 2
                vert = [rec.xmin(1) rec.xmin(2);
                        rec.xmin(1) rec.xmax(2);
                        rec.xmax(1) rec.xmax(2);
                        rec.xmax(1) rec.xmin(2)];
                return;
            end

            ind_fulldim = rec.getFullDims;
            n_fulldim = length(ind_fulldim);

            vert = repmat(rec.xmin, 2^n_fulldim, 1);

            for i = 2:2^n_fulldim
                act_fulldim = find(de2bi(i-1, n_fulldim));
                vert(i,ind_fulldim(act_fulldim)) = rec.xmax(ind_fulldim(act_fulldim));
            end
        end

        function [handle, HHH] = plotRec(rec, color, alpha, emph, multiple)
            if nargin<2
                color = [1 0 0];
            end
            if nargin<3
                alpha = 1;
            end
            if nargin<4
                emph = 0;
            end
            if nargin<5
                multiple = 1;
            end
            if length(rec)>1
                % We have an array
                hh = ishold; if (~hh) clf; end
                hold on
                if multiple
                    color = autumn(length(rec));
                end
                for i = 1:length(rec)
                    handle = plot(rec(i), color(min(i, size(color,1)),:), alpha);
                end
                if (~hh) hold off; end
                return;
            end	

            HHH = [];
            if length(rec) == 0
                disp('Warning: tried to plot nonexistent cell')
                % do nothing
                return;
            end

            if rec.dim>3
                error('Cant plot in dimension larger than 3')
            end

            if rec.isEmptySet
                return;
            end

            hh = ishold;
            if ~hh
                clf
            end
            hold on;

            if rec.dim==1
                h = line([rec.xmin rec.xmax], [0 0]);
                set(h, 'Color', color)
                HHH = [HHH,h];
            elseif rec.dim==2
                    vert = rec.getVertices;
                if rec.isFullDim
                    h = fill(vert(:,1), vert(:,2), 'r');
                    set(h, 'FaceColor', color)
                    set(h, 'FaceAlpha', alpha)
                    HHH = [HHH,h];
                else
                    h = line(vert(:,1), vert(:,2), 'r');
                    set(h, 'Color', color)
                    HHH = [HHH,h];
                end
            elseif rec.dim==3
                n_fulldim = length(rec.getFullDims);
                if n_fulldim == 3
                    projx = projection(rec, [2 3]);
                    vertx = projx.getVertices;
                    h = fill3(rec.xmin(1)*ones(4,1), vertx(:,1), vertx(:,2), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha); 
                    if (emph)
                        set(h,'Linewidth',1)
                    else
                        set(h,'EdgeColor','None')
                    end
                    HHH = [HHH,h];
                    h = fill3(rec.xmax(1)*ones(4,1), vertx(:,1), vertx(:,2), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha); 
                    if (emph)
                        set(h,'Linewidth',1)
                    else
                        set(h,'EdgeColor','None')
                    end
                    HHH = [HHH,h];
                    projy = projection(rec, [1 3]);
                    verty = projy.getVertices;
                    h = fill3(verty(:,1), rec.xmin(2)*ones(4,1), verty(:,2), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha); 
                    if (emph)
                        set(h,'Linewidth',1)
                    else
                        set(h,'EdgeColor','None')
                    end
                    HHH = [HHH,h];
                    h = fill3(verty(:,1), rec.xmax(2)*ones(4,1), verty(:,2), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha); 
                    if (emph)
                        set(h,'Linewidth',1)
                    else
                        set(h,'EdgeColor','None')
                    end
                    HHH = [HHH,h];
                    projz = projection(rec, [1 2]);
                    vertz = projz.getVertices;
                    h = fill3(vertz(:,1), vertz(:,2), rec.xmin(3)*ones(4,1), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha); 
                    if (emph)
                        set(h,'Linewidth',1)
                    else
                        set(h,'EdgeColor','None')
                    end
                    HHH = [HHH,h];
                    h = fill3(vertz(:,1), vertz(:,2), rec.xmax(3)*ones(4,1), 'r'); set(h, 'FaceColor', color); set(h, 'FaceAlpha', alpha); 
                    if (emph)
                        set(h,'Linewidth',1)
                    else
                        set(h,'EdgeColor','None')
                    end
                    HHH = [HHH,h];
                elseif n_fulldim == 2
                    ind_fulldims = rec.getFullDims;
                    ind_flatdim = rec.getFlatDims;
                    flatdimval = rec.xmin(ind_flatdim);
                    proj = projection(rec, ind_fulldims);
                    vert = getVertices(proj);
                    if ind_flatdim == 1
                        h = fill3(flatdimval*ones(4,1), vert(:,1), vert(:,2), 'r');
                    elseif ind_flatdim == 2
                        h = fill3(vert(:,1), flatdimval*ones(4,1), vert(:,2), 'r');
                    else
                        h = fill3(vert(:,1), vert(:,2), flatdimval*ones(4,1), 'r');
                    end
                    set(h, 'FaceColor', color)
                    set(h, 'FaceAlpha', alpha)
                    HHH = [HHH,h];
                else
                    vert = rec.getVertices;
                    h = line(vert(:,1), vert(:,2), vert(:,3));
                    set(h, 'Color', color)
                    HHH = [HHH,h];
                end
            end
            if ~hh
                hold off
            end
            handle = h;
        end
        
        function irec = intersect(rec1, rec2)
            % Computes the intersection between two Rec's rec1 and rec2

            if length(rec1)>1
                irec = [];
                for rec = rec1
                    isect = intersect(rec, rec2);
                    if ~isect.isEmptySet
                        irec = [irec isect];
                    end
                end
                return;
            end
            if length(rec2)>1
                irec = intersect(rec2, rec1);
                return;
            end

            imin = max(rec1.xmin, rec2.xmin);
            imax = min(rec1.xmax, rec2.xmax);

            % if ~isempty(setxor(rec1.ap, rec2.ap)) && (~isempty(rec1.ap) || ~isempty(rec2.ap))
            % 	warning('intersecting sets with different APs. Keeping both..')
            % end
            ap = union(rec1.ap, rec2.ap);
            irec = Rec([imin; imax], ap);
        end
        
        function irecs = multple_intersect_one(recs, rec)
            irecs = [];
            for i = 1:1:length(recs)
                irecs = [irecs, intersect(recs(i), rec)];
            end
        end
        
        function evolved_rec = evolve(rec, A,K)
            v_rec = my_exp(rec.xmin, rec.xmax);
            v_evolved_rec = A*v_rec + K;
            evolved_rec = Rec([min(v_evolved_rec,[],2)';...
                               max(v_evolved_rec,[],2)']);
        end
        
        function poly = Polyfy(rec)
            poly = [];
            for i = 1:1:length(rec)
                poly = [poly, Polyhedron('lb',rec(i).xmin,'ub',rec(i).xmax)];
            end
        end
	end
end


