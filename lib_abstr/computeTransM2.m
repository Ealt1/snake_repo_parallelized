function [cell_idx_all, trans, trans_out, trans_super] = computeTransM2...
    (act_set, gpart, new_cell_idx, cell_idx)
% Compute transition matrix for the partition part
% Returns a sparse matrix M where M(i,j)=1 if there
% is a transition from cell i to cell j.
% 
% Inputs: 
% Outputs:

N = length(new_cell_idx);
M = length(act_set);

if M>1
    % If there are several modes, 
    % compute transitions for each one
    trans = ones(N, N, M);
    trans_out = ones(N, 1, M);
    trans_super = ones(N, 1, M);
    for i=1:1:M
        i
        [~, trans(:,:,i), trans_out(:,:,i), trans_super(:,:,i)] = ...
            computeTransM2(act_set{i}, gpart, new_cell_idx, cell_idx);
    end
    cell_idx_all = new_cell_idx;
    return;
end


% trans = sparse(N, N);
% trans_out = sparse(N,1);
% trans_super = sparse(N,1);
trans = zeros(N, N);
trans_out = zeros(N,1);
trans_super = zeros(N,1);
for i = 1:1:N-1
    % i
    idx_i = new_cell_idx(i);
    adj_idx = gpart.getNeighbor_idx(idx_i);
    for idx_j=adj_idx
        if ~ismember(idx_j,cell_idx) 
            if gpart.isTrans(idx_i, idx_j, act_set)
                j = find(new_cell_idx == idx_j);
                trans(i,j) = 1;
                if isempty(j)
                    trans_out(i) = 1;
                end
            end
        else
            trans_super(i) = 1;
        end
    end
    
    if gpart.isTransOut(idx_i, act_set)
        trans_out(i) = 1;
    end
    
    % Hacking part
    % begin{hack}
    
    i_density = act_set.u(1);
    PB_heat = act_set.u(4);

    cord_i = gpart.idx2cord(idx_i);
    Rec_i = gpart.cord2Rec(cord_i);
    Ccell_Rec = Rec([Rec_i.xmin(1:7)',Rec_i.xmax(1:7)']);

    % 1) battery charge/doscharge requirement
    [mEFC_cell, ~] = optimumEFC_cell(act_set,Ccell_Rec);
    if 1.7e4 + PB_heat - i_density * mEFC_cell * 412 * 160 > 20*3.2^2/6*1e3
        trans_out(i) = 1;
    end
    
    % 1) FC output power must be positive
    if mEFC_cell < 0
        trans_out(i) = 1;
    end
    % end{hack}
    
    trans(i,i) = 1;
end
cell_idx_all = new_cell_idx;

