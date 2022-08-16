function frec = getFacetI(rec,i,j)
	% return facet number i of rec
    if (nargin < 3)
        j = [];
    end
	if abs(i)<1 || abs(i)>2*rec.dim
		error('index out of range')
	end
	fmin = rec.xmin;
	fmax = rec.xmax;
    
	if i < 0
		fmax(abs(i)) = fmin(abs(i));
	else
		fmin(abs(i)) = fmax(abs(i));
    end
    fmin(j) = [];
    fmax(j) = [];
	frec = Rec([fmin; fmax]);
end