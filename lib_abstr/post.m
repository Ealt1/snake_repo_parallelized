function post_set = post(states, transM, act)
  % POST: computes the 1-step Post set of a set of states under the action(s) act
  % if no action is specified, the Post set of all action is computed
	if nargin<3
		act = 1:size(transM,3);
	end
	post_set = [];
	for a = act
		post_set = union(post_set, find(sum(transM(states,:,a),1)>=1));
	end
	post_set = reshape(post_set, 1, length(post_set));
end
