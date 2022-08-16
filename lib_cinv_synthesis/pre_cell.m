function ret = pre_cell(states, trans_set, act)
	% PRE: computes the 1-step Pre set of a set of states under the action(s) act
	% if no action is specified, the Pre set of all action is computed
	if nargin<3
		act = 1:size(trans_set,3);
	end
	ret = [];
    A = zeros(act);
	for a = act
        if ~isempty(states)
            my_trans = sparse(trans_set(:,a,states));
            ret = union(ret, find(sum(my_trans,2)>=1));
        end
        %my_trans(:,a, states)
	end
	ret = reshape(ret, 1, length(ret));
    size(ret)
end