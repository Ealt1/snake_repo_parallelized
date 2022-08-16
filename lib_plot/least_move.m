function u_idx = least_move(x, uu_idx, gpart_u)

% ----------------- Heuristic 1 ----------------- 
move_size = [];
for u_idx = uu_idx
    u = gpart_u.idx2mid(u_idx);
    ms = 0; 
    if u(1) ~= x(1) 
        ms = ms + 1.5;
    end
    if u(2) ~= x(2) 
        ms = ms + 1;
    end
    if u(3) ~= x(3) 
        ms = ms + 1;
    end
    if ms == 0
        return;
    end
    move_size = [move_size, ms];
end
[~,i_least] = min(move_size); 
i_least = find(move_size == move_size(i_least));
if length(i_least) <= 1
    u_idx = uu_idx(i_least);
    return; 
end

charge_size = [];
for i = i_least
    u_idx = uu_idx(i);
    u = gpart_u.idx2mid(u_idx);
    cs = 0;
    if u(1) == u(2) 
        cs = cs + 1;
    end
    if u(1) == u(3) 
        cs = cs + 1;
    end
    charge_size = [charge_size, cs];
end
[~,i_least] = min(charge_size); 
u_idx = uu_idx(i_least);

% ----------------- Heuristic 2 ----------------- 

% u_idx_backup = uu_idx(1);
% flag = 0; 
% for u_idx = uu_idx
%     u = gpart_u.idx2mid(u_idx);
%     if ~flag
%         if u(2) == x(2) && u(3) == x(3)
%             u_idx_backup = u_idx;
%             flag = 1;
%         else
%             if u(2) == x(2) || u(3) == x(3)
%                 u_idx_backup = u_idx;
%             end
%         end
%     end
%     if u(1) == x(1)
%         if u(2) == x(2) && u(3) == x(3)
%             return;
%         else
%             if u(2) == x(2) || u(3) == x(3)
%                 return;
%             end
%         end
%     end
% end
% 
% u_idx = u_idx_backup;



