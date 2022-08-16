function X_safe = label_sys(sys)

X_safe = [];
for i = 1:1:length(sys.gpart_x)
    x = sys.gpart_x.idx2cord(i);
    if ismember(x(1), sys.p_large_safe)
        X_safe = [X_safe, i]; 
    end
end
