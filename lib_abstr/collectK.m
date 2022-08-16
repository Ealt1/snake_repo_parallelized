function k = collectK(cell_idx_idx, K)

list = {}
for i = cell_idx_idx
    if intersect(K{i},list)
end

