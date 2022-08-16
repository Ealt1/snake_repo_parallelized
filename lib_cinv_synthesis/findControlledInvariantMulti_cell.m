function [C,K] = findControlledInvariantMulti_cell(X, trans_set)

C = X;
while 1 % termination guaranteed - finite
	[Ct, K] = CPreMulti_cell(C, trans_set);
	if isempty(setdiff(C, Ct))
		% remove state/controller combinations which are not in C		
		% K(~ismember(Ct, C)) = [];
        K = K(ismember(Ct, C));
		Ct(~ismember(Ct, C)) = [];
		C = Ct;
		break;
	else
		C = intersect(Ct, C);
	end
end

if size(C,1)>1
    C = C';
end

