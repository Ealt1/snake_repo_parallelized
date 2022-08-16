function [cell_idx_all, trans, trans_out] = computeTransM_PowerCuitCostumized(act_set, gpart, cell_idx)
% Compute transition matrix for the partition part
% Returns a sparse matrix M where M(i,j)=1 if there
% is a transition from cell i to cell j.
% 
% Inputs: 
% Outputs:

% cell_idx_all = getNeighbor_idx(gpart, cell_idx, 1);
cell_idx_all = cell_idx;
N = length(cell_idx_all);
M = length(act_set);

if M>1
    % If there are several modes, 
    % compute transitions for each one
    trans = ones(N, N, M);
    trans_out = ones(N, 1, M);
    for i=1:1:M
        i
        [cell_idx_all, trans(:,:,i), trans_out(:,:,i)] = ...
            computeTransM_PowerCuitCostumized(act_set{i}, gpart, cell_idx);
    end
    return;
end


[single_A, single_K] = getAKd_PowerCircuit(1, 25e-9);
% trans = sparse(N, N);
% trans_out = sparse(N,1);
trans = zeros(N, N);
trans_out = zeros(N,1);
for i = 1:1:N
    % i
    cell_idx_i = cell_idx_all(i);
    cell_cord_i = gpart.idx2cord(cell_idx_i);

    % update: i,v, in rec
    cell_Rec_i = cord2Rec(gpart, cell_cord_i);
    cell_Rec_i_12 = Rec([cell_Rec_i.xmin(1:2)',cell_Rec_i.xmax(1:2)']'); 
    cell_Recs_j = evolve(cell_Rec_i_12, act_set.A, act_set.K);

    min_cell_cord = locatePoint_cord(gpart,[cell_Recs_j.xmin';1e-13;1.8]);
    max_cell_cord = locatePoint_cord(gpart,[cell_Recs_j.xmax';1e-13;1.8]);

    tro = ~contains(Rec([gpart.domain.xmin(1:2)',gpart.domain.xmax(1:2)']'), ...
                    Rec([cell_Recs_j.xmin(1:2)',cell_Recs_j.xmax(1:2)']'));

    nnn = gpart.dim;
    for iii = 1:1:nnn
        if min_cell_cord(iii) == -Inf
            min_cell_cord(iii) = 1;
        end
        if max_cell_cord(iii) == Inf
            max_cell_cord(iii) = length(gpart.grid{iii})-1;
        end
    end

    cell_cord_j_1 = min_cell_cord(1):1:max_cell_cord(1);
    cell_cord_j_2 = min_cell_cord(2):1:max_cell_cord(2);

    % determine which modes are available using the dewelling time status
    if cell_cord_i(4) == 1 % t = 0, can chose stay with mode 2, or switch to mode 1
        if ~any(act_set.K) % mode 2 (use the fact that K2 = [0;0]) 
            % update: t reset to 0, v^ fixed
            cell_cord_j_3 = cell_cord_i(3);
            cell_cord_j_4 = 1;
        else % mode 1
            % update: t
            cell_cord_j_4 = cell_cord_i(4) + act_set.t; 
            % update: v^
            if cell_cord_i(4) <= 4 && cell_cord_j_4 >= 5
                t_v_evolve2vhat = 5 - cell_cord_i(4);
                AA = eye(2);
                KK = [0;0];
                for tt = 1:1:t_v_evolve2vhat
                    AA = AA*single_A;
                    KK = single_A*KK + single_K;
                end
                % the old cell_Recs_j is useless at this point 
                cell_Recs_j = evolve(cell_Rec_i_12, AA, KK);

                % !!!!!!!!!!!!!!!!!!!!!!!!!
                % tro = tro || ~contains(Rec([gpart.domain.xmin(1:2)',gpart.domain.xmax(1:2)']'), ...
                %     Rec([cell_Recs_j.xmin(1:2)',cell_Recs_j.xmax(1:2)']'));
                    
                min_cell_cord = locatePoint_cord(gpart,[cell_Recs_j.xmin';1e-13;1.8]);
                max_cell_cord = locatePoint_cord(gpart,[cell_Recs_j.xmax';1e-13;1.8]);

                nnn = gpart.dim;
                for iii = 1:1:nnn
                    if min_cell_cord(iii) == -Inf
                        min_cell_cord(iii) = 1;
                    end
                    if max_cell_cord(iii) == Inf
                        max_cell_cord(iii) = length(gpart.grid{iii})-1;
                    end
                end
                
                cell_cord_j_3 = min_cell_cord(3):1:max_cell_cord(3);
            else
                cell_cord_j_3 = cell_cord_i(3); 
            end
        end
    else % t ~= 0, can only proceed with mode 1
        if ~any(act_set.K) % mode 2 (use the fact that K2 = [0;0]) 
            % mode 2 forbidden, -> out
            trans_out(i) = 1;
        else
            if cell_cord_i(4) + act_set.t > 11
                % forbidden, -> out
                trans_out(i) = 1;
            else % mode 1
                % update: t
                cell_cord_j_4 = cell_cord_i(4) + act_set.t; 
                % update: v^
                if cell_cord_i(4) <= 4 && cell_cord_j_4 >= 5
                    t_v_evolve2vhat = 5 - cell_cord_i(4);
                    AA = eye(2);
                    KK = [0;0];
                    for tt = 1:1:t_v_evolve2vhat
                        AA = AA*single_A;
                        KK = single_A*KK + single_K;
                    end
                    % the old cell_Recs_j is useless at this point 
                    cell_Recs_j = evolve(cell_Rec_i_12, AA, KK);

                    % !!!!!!!!!!!!!!!!!!!!!!!!!
                    % tro = tro || ~contains(Rec([gpart.domain.xmin(1:2)',gpart.domain.xmax(1:2)']'), ...
                    %                         Rec([cell_Recs_j.xmin(1:2)',cell_Recs_j.xmax(1:2)']'));

                    min_cell_cord = locatePoint_cord(gpart,[cell_Recs_j.xmin';1e-13;1.8]);
                    max_cell_cord = locatePoint_cord(gpart,[cell_Recs_j.xmax';1e-13;1.8]);

                    nnn = gpart.dim;
                    for iii = 1:1:nnn
                        if min_cell_cord(iii) == -Inf
                            min_cell_cord(iii) = 1;
                        end
                        if max_cell_cord(iii) == Inf
                            max_cell_cord(iii) = length(gpart.grid{iii})-1;
                        end
                    end

                    cell_cord_j_3 = min_cell_cord(3):1:max_cell_cord(3);
                else
                    cell_cord_j_3 = cell_cord_i(3); 
                end
            end
        end
    end

    
    % if gpart.isTransOut(cell_idx_i, act_set)
    %     trans_out(i) = 1;
    % end
    if tro
        trans_out(i) = 1;
    else
        % !!!!!!!!!!!!!!!!!!!!!!!!!!
        % min_cell_cords_j = [min(cell_cord_j_1); ...
        %                     min(cell_cord_j_2); ...
        %                     min(cell_cord_j_3); ...
        %                     min(cell_cord_j_4)];
        % max_cell_cords_j = [max(cell_cord_j_1); ...
        %                     max(cell_cord_j_2); ...
        %                     max(cell_cord_j_3); ...
        %                     max(cell_cord_j_4)];
        min_cell_cords_j = [min(cell_cord_j_1); ...
                            min(cell_cord_j_2); ...
                            cell_cord_i(3); ...
                            min(cell_cord_j_4)];
        max_cell_cords_j = [max(cell_cord_j_1); ...
                            max(cell_cord_j_2); ...
                            cell_cord_i(3); ...
                            max(cell_cord_j_4)];
        cell_cords_j = my_expNs(min_cell_cords_j, max_cell_cords_j, ...
            max_cell_cords_j - min_cell_cords_j + ones(length(min_cell_cords_j),1));
        cell_idxs_j = gpart.idx2cord(cell_cords_j);

        for idx_j=cell_idxs_j
            j = find(cell_idx_all == idx_j);
            trans(i,j) = 1;
            if isempty(j)
                trans_out(i) = 1;
            end
        end
    end
    
    % trans(i,i) = 1;
end

