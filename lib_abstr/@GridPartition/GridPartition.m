classdef GridPartition<handle

	properties (SetAccess=protected)
        domain;		% Domain (including ddim) 
		dim;		% Dimension of partition, (including ddim), 'n'
        grid        % Griding information (including ddim) 
	end

	methods
		function gpart = GridPartition(grid)
			% PARTITION: Create a GridPartition object
			%
			% SYNTAX
			% ------
			%
			%	gpart = GridPartition(grid)
			% 
			% INPUT
			% -----
			%	
            %   grid    1xn cell, e.g. {[1,2,4,5], ... ,[100, 200]}
			%
			% SEE ALSO
			% --------
			%            
            gpart.grid = grid;
            % sanity check
            if any(cellfun(@length,grid) <= 1)
                error('Invalid grid: unmber of grid points along each dimension must be >= 2')
            end
            
            gpart.dim = length(grid);
            
            xmin = cellfun(@min,grid); 
            xmax = cellfun(@max,grid); 
            gpart.domain = Rec([xmin; xmax]);
        end
        
        function l = length(obj)
        	% Overload length, i.e., number of cells          
            % length/resolution of grid along each dimension
            L = cellfun('length',obj.grid)-1; 
            l = prod(L);
        end
        
        function max_cord = get_max_cord(gpart, dim)
            L = cellfun('length',gpart.grid)-1; 
            max_cord = L(dim);
        end
        
        function min_mid_pnt = get_min_mid(gpart, dim)
            grid_dim = gpart.grid{dim};
            min_mid_pnt = (grid_dim(1) + grid_dim(2))/2;
        end

        function max_mid_pnt = get_max_mid(gpart, dim)
            max_cord = get_max_cord(gpart, dim);
            grid_dim = gpart.grid{dim};
            max_mid_pnt = (grid_dim(max_cord) + grid_dim(max_cord+1))/2;
        end
        
        function rv = isValid_cord(gpart, cell_cords)
            % return 0 if any of the coordinates are outside the domain,
            % and return 1 otherwise 
            m = size(cell_cords,2);
            L = cellfun('length',gpart.grid)-1;
            rv = ~(any(any(cell_cords >= repmat(L'+1,1,m))) ...
                   || any(any(cell_cords <= 0))); 
        end
        
        function ret = isContainedByDomain(gpart, rec)
            ret = contains(gpart.domain, rec);
        end
        
        function cell_idxs = cord2idx(gpart, cell_cords)
            % *
            % input: 
            %   cell_cords, nxm integer matrix
            % output:
            %   cell_idxs, 1xm integer row vector
            % assume:
            n = size(cell_cords,1);
            m = size(cell_cords,2);
            L = cellfun('length',gpart.grid)-1;
            W = [];
            for i = 1:1:n
                W = [W;prod(L(i+1:n))];
            end
            cei = ones(n,1);
            cei(n) = 0;
            cell_idxs = W'*(cell_cords - cei*ones(1,m));
        end
        
        function cell_cords = idx2cord(gpart, cell_idxs)
            % *
            % input: 
            %   cell_idxs, 1xm integer row vector
            % output:
            %   cell_cords, nxm integer matrix
            % assume:
            n = gpart.dim;
            m = size(cell_idxs,2);
            L = cellfun('length',gpart.grid)-1;
            cell_cords = [];
            for i = 1:1:n
                cell_cords = [cell_cords; ...
                    floor((cell_idxs-1)/prod(L(i+1:n)))+1];
                cell_idxs = cell_idxs - (cell_cords(end,:)-1)*prod(L(i+1:n));
            end
        end
        
        function cell_Recs = cord2Rec(gpart, cell_cords, poly_or_rec)
            % *
            % input: 
            %   cell_cords, nxm integer matrix
            % output:
            %   cell_Recs, 1xm Rec row vector
            % assume:
            
            if nargin < 3
                poly_or_rec = 'Rec';
            end
            cell_Recs = [];
            n = size(cell_cords,1);
            m = size(cell_cords,2);
            for j = 1:1:m
                xmin = [];
                xmax = [];
                for i = 1:1:n
                    xmin = [xmin; gpart.grid{i}(cell_cords(i,j))];
                    xmax = [xmax; gpart.grid{i}(cell_cords(i,j)+1)];
                end
                
                if strcmp(poly_or_rec, 'Rec')
                    cell_Recs = [cell_Recs, Rec([xmin';xmax'])];
                else
                    cell_Recs = [cell_Recs, Polyhedron('lb',xmin, 'ub', xmax)];
                end
            end
        end
        
        function cell_Recs = idx2Rec(gpart, cell_idxs, poly_or_rec)
            if nargin < 3
                poly_or_rec = 'Rec';
            end
            cell_cords = idx2cord(gpart, cell_idxs); 
            cell_Recs = cord2Rec(gpart, cell_cords, poly_or_rec); 
        end
        
        function cell_mid_pnts = idx2mid(gpart, cell_idxs)
            cell_Recs = gpart.idx2Rec(cell_idxs, 'Rec'); 
            cell_mid_pnts = getMidpoint(cell_Recs); 
        end
        
        function cell_idxs = pnt2idx(gpart, x)
            cell_idxs = locatePoint_idx(gpart,x); 
        end
        
        function Ncell_idxs = getFacetNeighbor_idx(gpart, cell_idxs)
            % return neigboring cells of 'cell_idxs', excluding everything 
            % in 'cell_idxs' 
            % input: 
            %   cell_idxs, 1xm integer row vector
            % output:
            %   Ncell_idxs, 1x? integer row vector
            % assume:
            % remark: 
            %   Since we are using rec-partition, non-facet adjacent 
            %   cases are not considered.
            
            cell_cords = gpart.idx2cord(cell_idxs); 
            Ncell_cords = [];
            
            n = gpart.dim;
            m = size(cell_cords,2);
            L = cellfun('length',gpart.grid)-1;
            
            for j = 1:1:m
                cell_cordj = cell_cords(:,j);
                for i = 1:1:n
                    cell_cordj_perturbi_plus = cell_cordj;
                    cell_cordj_perturbi_plus(i) = ...
                        cell_cordj_perturbi_plus(i) + 1;
                    cell_cordj_perturbi_minus = cell_cordj;
                    cell_cordj_perturbi_minus(i) = ...
                        cell_cordj_perturbi_minus(i) - 1;
                    if (cell_cordj_perturbi_plus(i) <= L(i))
                        Ncell_cords = [Ncell_cords, ... 
                        cell_cordj_perturbi_plus];
                    end
                    if (cell_cordj_perturbi_minus(i) >= 1)
                        Ncell_cords = [Ncell_cords, ...
                            cell_cordj_perturbi_minus];
                    end
                end
            end
            Ncell_idxs = gpart.cord2idx(Ncell_cords);
            Ncell_idxs = setdiff(Ncell_idxs,cell_idxs);
        end
        
        function Pcell_cords = perturbAll(gpart, cell_cord, Pplus, Pminus)
            % return "perturbed cells" in coordinates
            % used only by 'getNeighborAll'
            % input:
            %   cell_cord, nx1 integer col vector (only one cordinite)
            % output:
            %   Pcell_cord, nx? integer matrix
            % assume:
            %   cell_cord contains only one column.
            
            if (nargin <3)
                Pplus = ones(1,gpart.dim);
                Pminus = ones(1,gpart.dim);
            end
            
            nout = gpart.isNeighborOutside_cord(cell_cord);
            pplus = ~((nout>0)|(nout==2)); 
            pplus = pplus & Pplus;
            pminus = ~((nout<0)|(nout==2)); 
            pminus = pminus & Pminus;
            Pcell_cords = my_exp3...
                (cell_cord - pminus, cell_cord,cell_cord + pplus);           
            Pcell_idx = setdiff(gpart.cord2idx(Pcell_cords),...
                gpart.cord2idx(cell_cord));
            Pcell_cords = gpart.idx2cord(Pcell_idx);
        end 
        
        function NAcell_idxs = getNeighborAll_idx(gpart, cell_idxs, Pplus, Pminus)
            % return neigboring cells of 'cell_idx', excluding 'cell_idx' 
            % input: 
            %   cell_idx, 1xm integer row vector
            %   extend, bool scalar, 1: extend along dummy dimension, 0:
            %   otherwise
            % output:
            %   NAcell_idx, 1x? integer row vector
            % assume:
            % remark: 
            %   (1) different from getNeighbor in the following sense
            %   e.g., here [2 1 1 1] and [1 2 2 2] are considered to be 
            %   adjecent, but not in getNeighbor.
            if (nargin <3)
                Pplus = ones(1,gart.dim);
                Pminus = ones(1,gart.dim);
            end
            
            cell_cords = gpart.idx2cord(cell_idxs); 
            NAcell_cords = [];
            m = size(cell_cords,2);
            for j = 1:1:m
                cell_cordj = cell_cords(:,j);
                NAcell_cords = [NAcell_cords, ... 
                gpart.perturbAll(cell_cordj, Pplus, Pminus)];
            end
            NAcell_idxs = gpart.cord2idx(NAcell_cords);
            NAcell_idxs = setdiff(NAcell_idxs,cell_idxs);
        end
        
        function Ncell_idx = getNeighborWithin_idx...
                (gpart, cell_idx, cellwithin_idx)
            % ?
            % return neigboring cells of 'cell_idx', excluding cells 
            % outside 'cellwithin_idx' 
            % used only by 'isBoundaryCellWithin'
            Ncell_idx = gpart.getFacetNeighbor_idx(cell_idx);
            Ncell_idx = intersect(Ncell_idx, cellwithin_idx);
        end
        
        function bc = isBoundaryCellWithin...
                (gpart, cell_idx, cellwithin_idx)
            % input:
            %   cell_idx, 1xm integer row vector
            % output:
            %   bc, 1xm bool row vector, 1: cellI
            L = cellfun('length',gpart.grid)-1;
            bc = ones(1,length(cell_idx));
            for i = 1:1:length(cell_idx)
                l1 = length(gpart.getNeighborWithin_idx...
                    (cell_idx(i), cellwithin_idx, extend));
                L_idx = 1:1:length(L);
                l2 = length(find(L(L_idx)==1));

                if l1 == 2*(length(L_idx)-l2)
                    bc(i) = 0;
                end
            end
        end 
        
        function nout = isNeighborOutside_cord(gpart, cell_cords)
            % *
            % input: 
            %   cell_cords, nxm integer col vector
            % output:
            %   nout, nxm integer row vector, adjacent dimension,
            %         e.g., 1: c1 < out
            %              -1: c1 > out
            %               2: c1 >< out
            %               0: not adjacent to outside
            % assume: inputs contains only one coordinate
            
            if ~gpart.isValid_cord(cell_cords)
                error('Invalid coordinates');
            end
            n = gpart.dim;
            m = size(cell_cords,2);
            L = cellfun('length',gpart.grid)-1;
            nout_p1 = 3*(cell_cords == L'*ones(1,m));
            nout_m1 = -(cell_cords == ones(n,m));
            nout = nout_p1 + nout_m1;
            nout(find(nout == 3)) = 1;
        end
        
        function n = isNeighbor_cord(gpart, cell1_cord, cell2_cord)
            % input: 
            %   cell1,2_cord, nx1 integer col vector
            % output:
            %   n, boolean scalar, 1 for neighboring, 0 o.w.
            % assume: 
            %   c1, c2 are different
            if ~gpart.isValid_cord(cell1_cord) || ~gpart.isValid_cord(cell2_cord)
                error('Invalid coordinates');
            end
            if all(cell1_cord - cell2_cord < 1e-5)
                error('Two cells need to be different');
            end
            rdim = 1:1:gpart.dim;
            n = (norm(cell1_cord(rdim) - cell2_cord(rdim),1)<=1);
        end
        
        function adim = getAdjacentDim_cord...
                (gpart, cell1_cord, cell2_cord)
            % input: 
            %   cell1,2_cord, nx1 integer col vectors
            %   order of cell1,2 matters
            % output:
            %   adim, integer, adjacent dimension, 
            %         adim>0: normal c1->c2 is pointing positive direction
            %         e.g., 5: c1 < c2 
            %              -5: c1 > c2
            % assume: 
            %   c1 c2 are indeed adjacent
            
            if ~isNeighbor_cord(gpart, cell1_cord, cell2_cord)
                error('Two cells need to be adjacent');
            end
            % get "extended" cell_Recs
            cell1_Rec = cord2Rec(gpart, cell1_cord);
            cell2_Rec = cord2Rec(gpart, cell2_cord);
            I_Rec = intersect(cell1_Rec, cell2_Rec);
            adim = getFlatDims(I_Rec);
            % adim = setdiff(adim, gpart.ddim);
            adim = adim*sign(...
                ((cell2_Rec.xmax(adim) - cell1_Rec.xmax(adim))>=0)-0.5...
                );
        end
        
        function cell_2trans2_idxs = TransTo(gpart, cell_idx, fx)
            % return true if there is a flow from c1->c2 under fx
            % inputs: 
            %   fx, SASys
            % output:
            %   tr, boolean scalar
            cell_cord = gpart.idx2cord(cell_idx);
            cell_Rec = cord2Rec(gpart, cell_cord);
            evolved_cell_Rec = evolve(cell_Rec, fx.A, fx.K);
            cell_2trans2_idxs = getBox_idx(gpart,evolved_cell_Rec);
        end
        
        function cell_2trans2_idxs = TransToBisim(gpart, cell_idx, fx)
            % return true if there is a flow from c1->c2 under fx
            % inputs: 
            %   fx, SASys
            % output:
            %   tr, boolean scalar
            cell_cord = gpart.idx2cord(cell_idx);
            cell_Rec = cord2Rec(gpart, cell_cord);
            % evolved_cell_Rec = evolve(cell_Rec, fx.A, fx.K);
            x_center = (cell_Rec.xmin + cell_Rec.xmax)/2;
            xx = fx.A*x_center' + fx.K; 
            cell_2trans2_idxs = gpart.cord2idx(locatePoint_cord(gpart,xx));
        end
        
        function tr = isTrans(gpart, cell1_idx, cell2_idx, fx)
            % return true if there is a flow from c1->c2 under fx
            % inputs: 
            %   fx, SASys
            % output:
            %   tr, boolean scalar
            cell1_cord = gpart.idx2cord(cell1_idx);
            cell2_cord = gpart.idx2cord(cell2_idx);
            
            
            cell1_Rec = cord2Rec(gpart, cell1_cord);
            evolved_cell1_Rec = evolve(cell1_Rec, fx.A, fx.K);
            cell2_Rec = cord2Rec(gpart, cell2_cord);
            tr = isEmptySet(intersect(evolved_cell1_Rec, cell2_Rec));
            
            % adim = getAdjacentDim_cord(gpart, cell1_cord, cell2_cord);
            % adjFacet = getFacetI(cell2_Rec,-adim);
            % [m, M] = optimum_f(fx,adjFacet,adim);
            % tr = (M>=1e-10);
        end
        
        function tro = isTransOut(gpart, cell_idx, fx)
           % Inputs: fx, FCTMSys
           % Assume: cell_idx is scalar
           cell_cord = gpart.idx2cord(cell_idx);
           cell_Rec = cord2Rec(gpart, cell_cord);
           evolved_cell_Rec = evolve(cell_Rec, fx.A, fx.K);
           tro = ~contains(gpart.domain, evolved_cell_Rec);
        end
        
        function cell_idx = locatePoint_idx(gpart,x)
            % assume:
            %   state 'x' is not at boundary between two cells.
            n = gpart.dim;
            cell_cord = [];
            for i = 1:1:n
                cord_i = find(x(i) <= gpart.grid{i});
                % x(i)
                % gpart.grid{i}
                % cord_i = find((gpart.grid{i} <= x(i)));
                % if ((~isempty(cord_i))...
                %         &&(gpart.grid{i}(end) > x(i)))
                if ((~isempty(cord_i)) && (gpart.grid{i}(1)<=x(i)))
                    % cord_i(1)-1
                    cell_cord = [cell_cord; cord_i(1)-1];
                else
                    cell_idx = [];
                    return;
                end
            end
            cell_idx = gpart.cord2idx(cell_cord);
        end
        
        function cell_cord = locatePoint_cord(gpart,x)
            % assume:
            %   state 'x' is not at boundary between two cells.
            n = gpart.dim;
            cell_cord = [];
            for i = 1:1:n
                cord_i = find(x(i) <= gpart.grid{i});
                if x(i) <= gpart.grid{i}(1)
                % if x(i) < gpart.grid{i}(1) % for some reason, leads to
                % negative idx
                    cell_cord = [cell_cord; -Inf];
                else
                    if x(i) >= gpart.grid{i}(end)
                    % if x(i) > gpart.grid{i}(end)
                        cell_cord = [cell_cord; +Inf];
                    else
                        cord_i = find(x(i) <= gpart.grid{i});
                        cell_cord = [cell_cord; cord_i(1)-1];
                    end
                end
            end
        end
        
        function cell_idxs = locateBox_idx(gpart,xBox)
            % input: xBox, Rec, 
            % output:
            % assume: 
            min_cell_idx = locateState_idx(gpart,xBox.xmin');
            min_cell_cord = gpart.idx2cord(min_cell_idx);
            max_cell_idx = locateState_idx(gpart,xBox.xmax');
            max_cell_cord = gpart.idx2cord(max_cell_idx);
            
            % max_cell_cord - min_cell_cord + ones(length(min_cell_cord),1)
            cell_cords = my_expNs(min_cell_cord, max_cell_cord, ...
                max_cell_cord - min_cell_cord + ones(length(min_cell_cord),1));
            
            cell_idxs = gpart.cord2idx(cell_cords);
            cell_idxs = setdiff(cell_idxs,[]);
        end
        
        function cell_idxs = getBox_idx(gpart,xBox)
            % input: xBox, Rec, 
            % output:
            % assume: 
            min_cell_cord = locatePoint_cord(gpart,xBox.xmin');
            max_cell_cord = locatePoint_cord(gpart,xBox.xmax');
            if ismember(Inf, min_cell_cord)||ismember(-Inf, max_cell_cord)
                cell_idxs = [];
                return;
            end
            n = gpart.dim;
            for i = 1:1:n
                if min_cell_cord(i) == -Inf
                    min_cell_cord(i) = 1;
                end
                if max_cell_cord(i) == Inf
                    max_cell_cord(i) = length(gpart.grid{i})-1;
                end
            end
            cell_cords = my_expNs(min_cell_cord, max_cell_cord, ...
                max_cell_cord - min_cell_cord + ones(length(min_cell_cord),1));
            
            cell_idxs = gpart.cord2idx(cell_cords);
            cell_idxs = setdiff(cell_idxs,[]);
        end
   
        function plotCell_by_idx(gpart, cell_idxs, dx, dy, clr)
            % input: 
            %   cell_idxs, 1xm integer row vector
            %   dx, dy: projection dimensions
            if nargin < 5
                clr = [0;0;0];
            end
            
            cell_cords = gpart.idx2cord(cell_idxs);
            cell_Recs = gpart.cord2Rec(cell_cords);
            
            for i = 1:1:length(cell_Recs)
                Ri = cell_Recs(i);
                xmin = Ri.xmin(dx);
                xmax = Ri.xmax(dx);
                ymin = Ri.xmin(dy);
                ymax = Ri.xmax(dy);
                hold on;
                plot([xmin,xmax,xmax,xmin,xmin],...
                     [ymin,ymin,ymax,ymax,ymin],'Linewidth',1,'Color', clr,'linestyle','--'); 
            end
        end
        
        function plotGrid(gpart, dx, dy)
            % input: 
            %   dx, dy: projection dimensions
            
            X = gpart.grid{dx};
            Y = gpart.grid{dy};
            for x = X
                plot([x,x],[Y(1),Y(end)],'k');
                hold on
            end
            for y = Y
                plot([X(1),X(end)],[y,y],'k');
                hold on
            end
            axis([X(1)-0.2*(X(end)-X(1)), X(end)+0.2*(X(end)-X(1)),...
                  Y(1)-0.2*(Y(end)-Y(1)), Y(end)+0.2*(Y(end)-Y(1))]);
            % axis square;
        end
        
        function plotMidPnts(gpart, cell_idxs, lw, clr)
            if gpart.dim > 3
                error('dim > 3'); 
            end
            if nargin < 4
                clr = [0 0 0]; 
            end
            if nargin < 3
                lw = 1; 
            end
            mid = gpart.idx2mid(cell_idxs); 
            if gpart.dim == 3
                scatter3(mid(1,:), mid(2,:), mid(3,:),3,clr, 'Linewidth', lw); 
            else
                scatter(mid(1,:), mid(2,:), 3, clr, 'Linewidth', lw); 
            end
            hold on; 
        end
        
        function plotArrow(gpart, cell1_idx, cell2_idx, end_txt, clr)
            if nargin < 5
                clr = 'black';
            end
            if nargin < 4
                end_txt = '';
            end
            x1 = gpart.idx2mid(cell1_idx); 
            x2 = gpart.idx2mid(cell2_idx); 
            plot([x1(1),x2(1)], [x1(2),x2(2)], 'Linewidth', 1,'Color',clr,'Linestyle',':');
            hold on; 
            a = [cos(pi/6), -sin(pi/6); sin(pi/6), cos(pi/6)]*(x1-x2)*0.3; 
            b = [cos(-pi/6), -sin(-pi/6); sin(-pi/6), cos(-pi/6)]*(x1-x2)*0.3; 
            % plot([x2(1),x2(1)+a(1)], [x2(2),x2(2)+a(2)], 'Linewidth', 1,'Color',clr);
            % plot([x2(1),x2(1)+b(1)], [x2(2),x2(2)+b(2)], 'Linewidth', 1,'Color',clr);
            text(x2(1)-0.3,x2(2)-0.5,end_txt,'Interpreter','latex');
        end
        
        function addTxt(gpart, cell_idx, txt, dir, clr, fnt_sz)
            if nargin < 6
                fnt_sz = 15;
            end
            if nargin < 5
                clr = [0,0,0];
            end
            x = gpart.idx2mid(cell_idx); 
            hold on; 
            if dir == 'm'
                text(x(1),x(2)-0.2,txt,'Fontsize',fnt_sz, 'Color', clr);
            end
            if dir == 's'
                text(x(1)-0.3,x(2)-0.5,txt,'Interpreter','latex','Fontsize',fnt_sz, 'Color', clr);
            end
            if dir == 'n'
                text(x(1)-0.3,x(2)+0.5,txt,'Interpreter','latex','Fontsize',fnt_sz, 'Color', clr);
            end
            if dir == 'w'
                text(x(1)-1,x(2),txt,'Interpreter','latex','Fontsize',fnt_sz, 'Color', clr);
            end
            if dir == 'e'
                text(x(1)+0.3,x(2),txt,'Interpreter','latex','Fontsize',fnt_sz, 'Color', clr);
            end
        end
    end
end


