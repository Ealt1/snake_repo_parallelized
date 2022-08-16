function pre_set = pre_action_forall_cell(X, trans_set, act)
	% Computes the set pre_set, such that for all s \in pre_set,
	% s -> X under action act.
	N_states = size(trans_set(:,:,1),1);
	pre_cand = pre_cell(X, trans_set, act);
	bad_states = setdiff(post_cell(pre_cand, trans_set, act), X);
    if ~isempty(pre_cand) && ~isempty(bad_states)
        my_trans = sparse(trans_set(pre_cand,act,bad_states));
        pre_set = pre_cand(find(sum(my_trans),2) == 0);
    else
        pre_set = [];
    end
	pre_set = reshape(pre_set, 1, length(pre_set));
end