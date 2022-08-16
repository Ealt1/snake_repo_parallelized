classdef GridSample<handle

	properties (SetAccess=protected)
        domain;		 % Domain (including ddim) 
		dim;		 % Dimension of partition, (including ddim), 'n'
        % ddim % Dimension of discrete control
        grid;        % Griding information (including ddim) 
        pnt_list;    % list of actions (fx's), length 'M'
	end

	methods
		function gsamp = GridSample(domain, grid, ddim)
			% PARTITION: Create a GridPartition object
			%
			% SYNTAX
			% ------
			%
			%	gpart = GridPartition(domain, grid, ddim)
			% 
			% INPUT
			% -----
			%	
			%	domain 	Rec, domain of the Partition
            %   grid    1xn cell, e.g. {[1,2,4,5], ... ,[100, 200]}
            %   assume  domain.dim = length(grid) = n
            %   ddim    1x? vector, subset of {1,2,...,n}
			%
			% SEE ALSO
			% --------
			%
            gsamp.dim = domain.dim;
            % gsamp.ddim = ddim;
            gsamp.grid = grid;
            gsamp.domain = domain;
            
            % length/resolution of grid along each dimension
            L = cellfun('length',gsamp.grid);
            
            M = prod(L);
            gsamp.pnt_list = 1:M; 
        end
        
        function pnt_idx = cord2idx(gsamp, pnt_cord)
            % input: 
            %   cell_cord, nxm integer matrix
            % output:
            %   cell_idx, 1xm integer row vector
            % assume:
            n = size(pnt_cord,1);
            m = size(pnt_cord,2);
            L = cellfun('length',gsamp.grid);
            W = [];
            for i = 1:1:n
                W = [W;prod(L(i+1:n))];
            end
            cei = ones(n,1);
            cei(n) = 0;
            pnt_idx = W'*(pnt_cord - cei*ones(1,m));
        end
        
        function pnt_cord = idx2cord(gsamp, pnt_idx)
            % input: 
            %   cell_idx, 1xm integer row vector
            % output:
            %   cell_cord, nxm integer matrix
            % assume:
            n = gsamp.dim;
            m = size(pnt_idx,2);
            L = cellfun('length',gsamp.grid);
            pnt_cord = [];
            for i = 1:1:n
                pnt_cord = [pnt_cord; ...
                    floor((pnt_idx-1)/prod(L(i+1:n)))+1];
                pnt_idx = pnt_idx - (pnt_cord(end,:)-1)*prod(L(i+1:n));
            end
        end
        
        function [low_pnt, up_pnt] = findLargestRec(u_gsamp, K_idx)
            M = length(K_idx);
            K_cord = u_gsamp.idx2cord(K_idx);
            LU_cord_cand = [];
            count_old = 0;
            for ii = 1:1:M
                for jj = ii:1:M
                    cord_ii = K_cord(:,ii); 
                    cord_jj = K_cord(:,jj); 
                    L_cord_ij = min([cord_ii,cord_jj],[],2);
                    U_cord_ij = max([cord_ii,cord_jj],[],2);
                    if prod(U_cord_ij+1 - L_cord_ij) > count_old
                        count = 0;
                        for kk = 1:1:M
                            if all( L_cord_ij <= K_cord(:,kk) ...
                                    & K_cord(:,kk) <= U_cord_ij)
                                count = count + 1;
                            end
                        end
                        if count == prod(U_cord_ij+1 - L_cord_ij) % #pnts in R_ij
                            if count > count_old
                                LU_cord_cand = [L_cord_ij, U_cord_ij];
                                count_old = count; 
                            end
                        end
                    end
                end
            end
            [low_pnt, up_pnt] = ...
                u_gsamp.getRec(LU_cord_cand(:,1), LU_cord_cand(:,2));
        end

        function [low_pnt, up_pnt] = getRec(gsamp, low_cord, up_cord)
            low_pnt = [];
            up_pnt = [];
            for d = 1:1:length(gsamp.grid)
                low_pnt = [low_pnt; gsamp.grid{d}(low_cord(d))];
                up_pnt = [up_pnt; gsamp.grid{d}(up_cord(d))];
            end
        end
    end
end


