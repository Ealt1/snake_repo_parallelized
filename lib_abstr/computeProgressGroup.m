function calG = computeProgressGroup(abstr, cell_idx, domain_cell_idx)
    m = length(cell_idx);
    calG = cell(1,length(abstr.act_set));
    for k = 1:length(abstr.act_set)
        k
        calG{k} = {};
        for j = 1:m
            idx_j = cell_idx(j);
            if ~ismember(idx_j, setdiff(cell2mat(calG{k}),[]))
                [ts, ~] = abstr.gpart.isTransient(idx_j, abstr.act_set{k});
                if ts
                    [gcell_idx, ~] = expandProgressGroup(abstr.gpart, ...
                        idx_j, abstr.act_set{k}, domain_cell_idx);
                    calG{k} = [calG{k}, [gcell_idx]];
                end
            end
        end
    end
end