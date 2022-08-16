function trans_set = computeTransM_linear(act_set, gpart)
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
        trans_set_i = computeTransM_linear(act_set{i}, gpart);
        trans_set{end+1} = trans_set_i{1}; 
    end
    return;
end


% trans = sparse(N, N);
% trans_out = sparse(N,1);
trans_set = sparse(N+1, N+1);
for i = 1:1:N
    % i
    cell_idx_i = cell_idx(i);
    cell_cord_i = gpart.idx2cord(cell_idx_i);
    cell_Rec_i = cord2Rec(gpart, cell_cord_i, 'Polyhedron');
    % maybe implemented by affineMap 
    Post_Rec_i = cell_Rec_i.affineMap(act_set.A) + act_set.K + act_set.W.affineMap(act_set.E); 
    % Post_Rec_i = Polyhedron('V',(act_set.A*cell_Rec_i.V' + act_set.K*ones(1,size(cell_Rec_i.V,1)))') + Polyhedron('V',(act_set.E*act_set.W.V')');  
    Post_Rec_i_xmin = min(Post_Rec_i.V,[],1)';
    Post_Rec_i_xmax = max(Post_Rec_i.V,[],1)';
    
    tro = ~contains(Rec([gpart.domain.xmin',gpart.domain.xmax']), ...
                    Rec([Post_Rec_i_xmin,Post_Rec_i_xmax]));
    
    if tro
        trans_set(i,N+1) = 1;
    else 
        cell_idx_j = gpart.getBox_idx(Rec([Post_Rec_i_xmin,Post_Rec_i_xmax]));
        
        % 1) less conservative, but may experience cddmex error
        % cell_cord_j = gpart.idx2cord(cell_idx_j);
        % cell_Rec_j = gpart.cord2Rec(cell_cord_j, 'Polyhedron'); 
        % for jj = 1:1:length(cell_idx_j)
        %     j = cell_idx_j(jj); 
        %     if i == 1289 && j == 138
        %         xxx = 999;
        %     end
        %     if ~isEmptySet(intersect(cell_Rec_j(jj), Post_Rec_i))
        %         trans_set(i,j) = 1;
        %     end
        % end
        
        % 2) more conservative, not MPT, faster
        trans_set(i,cell_idx_j) = 1;
    end
end

trans_set(N+1, N+1) = 1;

trans_set = {trans_set}; 




