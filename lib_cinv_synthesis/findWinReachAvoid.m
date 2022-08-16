function BR = findWinReachAvoid(X_target, X_unsafe, trans_set)

BR1 = X_target; 
BR0 = [];
while length(BR0) < length(BR1) 
    BR0 = BR1; 
    [BR_new, ~] = CPreMulti_cell(BR0, trans_set);
    BR1 = setdiff(union(BR0, BR_new), X_unsafe);
end
BR = BR1; 


