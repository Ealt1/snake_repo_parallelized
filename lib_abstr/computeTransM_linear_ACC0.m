function trans_set = computeTransM_linear_ACC(act_set, gpart)
% Compute transition matrix for the partition part
% Returns a sparse matrix M where M(i,j)=1 if there
% is a transition from cell i to cell j.
% 
% Inputs: 
%   act_set, cell array of acts (A,K,E,W)
% Outputs:

cell_idx = 1:1:length(gpart); 
N = length(cell_idx);
M = length(act_set);

if M>1
    % If there are several modes, 
    % compute transitions for each one
    trans_set = {};
    for i=1:1:M
        i
        trans_set_i = computeTransM_linear_ACC(act_set(i), gpart);
        trans_set{end+1} = trans_set_i{1}; 
    end
    return;
else
    act_set = act_set{1}; 
end

A = act_set.A; 
K = act_set.K;
E = act_set.E; 
W = act_set.W; 

X_safe = Polyhedron('A',[1,-1/1.7, 0; [1 0 0]; -[0 1 0]],'b', [0; [25]; -[4]]); 
% X_safe = Polyhedron('A',[1,-1/1.7, 0; [1 0 0; 0 0 1]; -[0 1 0]],'b', [0; [25; 25]; -[4]]); 
X_safe_pseudo = Polyhedron('A',[0 0 -1],'b', [-1e-6]); 

trans_set = sparse(N+1, N+1);
for i = 1:1:N
    % i
    
    Rec_i = gpart.idx2Rec(i, 'Polyhedron'); 
    if ~contains(X_safe, Rec_i)
        tro = 1; 
    else
        post_Poly_i = Rec_i.affineMap(A) + K + W.affineMap(E);
        % mid_pnt_i = gpart.idx2mid(i);
        % post_Poly_i = Polyhedron('lb', A*mid_pnt_i + K - 1e-6*ones(3,1), 'ub', A*mid_pnt_i + K + 1e-6*ones(3,1)); 
        
        post_Poly_i = intersect(X_safe_pseudo, post_Poly_i); 
        % post_Poly_i = max(post_Poly_i, [1e-6; -Inf; 1e-6]);
        % post_Poly_i = min(post_Poly_i, [Inf; 100-1e-6; Inf]);

        if  (isEmptySet(post_Poly_i))
            xxx = 999; 
        end
        tro = (~contains(X_safe, post_Poly_i));
    end
    
    if tro
        trans_set(i,N+1) = 1;
    else
        post_Rec_i_xmin = min(post_Poly_i.V,[],1)';
        post_Rec_i_xmax = max(post_Poly_i.V,[],1)';
        cell_idx_j = gpart.getBox_idx(Rec([post_Rec_i_xmin,post_Rec_i_xmax])); 
        if isempty(cell_idx_j)
            trans_set(i,N+1) = 1;
        end
        for jj = 1:1:length(cell_idx_j)
            j = cell_idx_j(jj); 
            Rec_j = gpart.idx2Rec(j, 'Polyhedron'); 
            if ~isEmptySet(intersect(Rec_j, post_Poly_i))
                trans_set(i,j) = 1;
            end
        end
        
        % cell_idx_j = gpart.locatePoint_idx(post_Poly_i); 
        % trans_set(i,cell_idx_j) = 1;
    end
end

trans_set(N+1, N+1) = 1;

trans_set = {trans_set}; 




