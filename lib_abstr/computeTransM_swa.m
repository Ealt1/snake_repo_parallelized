function trans_set = computeTransM_swa(sys, gpart)
% Compute transition matrix for the partition part
% Returns a sparse matrix M where M(i,j)=1 if there
% is a transition from cell i to cell j.
% 
% Inputs: 
% Outputs:

cell_idx = 1:1:length(gpart); 

N = length(cell_idx); 
M = length(sys.U);
nx = size(sys.A{1},1);

trans_set = {};
for k=1:1:M
    k
    trans_set_k = sparse(N+1, N+1);
    Wk = Polyhedron('A',sys.XW{k}.A(:,nx+1:end),'b',sys.XW{k}.b);
    for i = 1:1:N
        % i
        cell_idx_i = cell_idx(i);
        cell_cord_i = gpart.idx2cord(cell_idx_i);
        cell_Rec_i = cord2Rec(gpart, cell_cord_i, 'Polyhedron');
        % maybe implemented by affineMap 
        
        Post_Poly_i = cell_Rec_i.affineMap(sys.A{k}) + sys.K{k} + Wk.affineMap(sys.E{k});  
        Post_Rec_i_xmin = min(Post_Poly_i.V,[],1)';
        Post_Rec_i_xmax = max(Post_Poly_i.V,[],1)';

        Post_Pec_i = Rec([Post_Rec_i_xmin'; Post_Rec_i_xmax']); 
        tro = ~contains(gpart.domain,Post_Pec_i);
        if tro
            trans_set_k(i,N+1) = 1;
        else 
            % 2D case
            if isEmptySet(Post_Pec_i)
                error('empty Rec j');
            end
            cell_idx_j = gpart.getBox_idx(Post_Pec_i);
            for j = cell_idx_j
                if ~isEmptySet(intersect(Post_Poly_i,cord2Rec(gpart, gpart.idx2cord(j), 'Polyhedron')))
                    trans_set_k(i,j) = 1;
                end
            end
        end
    end

    trans_set_k(N+1, N+1) = 1;
    trans_set{end+1} = trans_set_k; 
end




