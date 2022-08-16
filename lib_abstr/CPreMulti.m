function [pre_set, K] = CPreMulti(X, trans_set)
pre_set = [];
K = {};
for a = 1:size(trans_set, 3)
   pre_set_a = pre_action_forall(X, trans_set, a);
   new_states = setdiff(pre_set_a, pre_set);

   if ~isempty(intersect(pre_set_a,pre_set))
       for state = intersect(pre_set_a,pre_set)
          K{find(pre_set == state)} = [K{find(pre_set == state)}, a];
       end
   end
   for nstate = new_states
      K = {K{1:end}, a};
   end
   pre_set = [pre_set, new_states];

end
pre_set = reshape(pre_set, 1, length(pre_set));