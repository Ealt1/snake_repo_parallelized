classdef Abstraction<handle

	properties (SetAccess=protected)
		gpart;		% Grid-based partition
        act_set;    % Set of ations
        label;      % Labeling function
	end

	methods
		function abstr = Abstraction(gpart,act_set)
			% PARTITION: Create a Partition object
            % Remark: this is not a general tool, it is specitially
            %   desinged for FCTM problem.
			%
			% SYNTAX
			% ------
			%
			%	abstr = Abstraction(gpart,act_set)
			% 
			% INPUT
			% -----
			%	
			%   gpart, GridPartition
            %   act_set, cells of FCTMSys
            %   label, cells of aps, 
            %          length(label) = length(gpart.cell_list).
            %
			% SEE ALSO
			% --------
			%
            abstr.gpart = gpart;
            abstr.act_set = act_set;
            abstr.label = cell(1,length(gpart.cell_list));
        end

        function generateLabeling(abstr, everything_is_safe)
            % generate abstr.label
            % ap 
            %  1: goal
            %  2: unsafe
            %  5,6,7: reach_1,2,3
            if nargin < 2
                everything_is_safe = 0;
            end
            if everything_is_safe
                for cell_idx = 1:1:length(abstr.gpart.cell_list)
                    abstr.label{cell_idx} = union(...
                        abstr.label{cell_idx}, 1);
                end
            else
                for cell_idx = 1:1:length(abstr.gpart.cell_list)
                    % cell_idx
                    cell_Rec = cord2Rec(abstr.gpart, ...
                        idx2cord(abstr.gpart, cell_idx));
                    if (goal(cell_Rec))
                        abstr.label{cell_idx} = union(...
                            abstr.label{cell_idx}, 1);
                    end
                end
            end
        end
        
        function TransM = generateTransM(abstr)
            % build transM within cell_idx, the only other node is "out"
            % input: 
            % output: 
            % assume:
            cell_idx = abstr.gpart.cell_list;
            [~, TransM, trans_out] = computeTransM_PowerCuitCostumizedxxx...
                (abstr.act_set, abstr.gpart, cell_idx);
            TransM = [TransM, trans_out];
        end
        
        function calG = generateProgressGroup(abstr, cell_idx)
            % inputs:
            % outputs: 
            calG = cell(1,length(abstr.act_set));
            for k = 1:length(abstr.act_set)
                k
                Gk = {};
                for idx = cell_idx
                    if  ~ismember(idx,cell2mat(Gk))
                        [g, nvdim] = expandProgressGroup...
                            (abstr.gpart, idx, abstr.act_set{k});
                        if ~isempty(g)
                            Gk{end+1} = g;
                        end
                    end
                end
                calG{k} = Gk;
            end
        end
        
        function cell_idx_ap = getap(abstr, ap, cell_idx_domain)
            % ap scalar
            if (nargin < 3)
                cell_idx_domain = abstr.gpart.cell_list;
            end
            Loc = cellfun(@(x) x==ap, abstr.label,'Un',0);
            Loc = cellfun(@(x) any(x(:)),Loc);
            cell_idx_ap = find(Loc);
            cell_idx_ap = intersect(cell_idx_ap, cell_idx_domain);
        end
    end
end