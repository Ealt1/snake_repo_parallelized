function sys = load_sys6(flag)

if nargin < 1
    flag = '';
end

n = 7;
Map = zeros(n,n);
Map(1,4) = 1;
Map(4,5) = 1;
Map(4,6) = 1;
Map(5,2) = 1;
Map(5,6) = 1;
Map(6,7) = 1;
Map(7,3) = 1;
Map = (Map + Map' + eye(n)) > 0.5; 

sys.M = Map;

sys.p_large_safe = [4 5 6 7];
sys.p_target = [1 2 3];
l = length(sys.p_target);

sys.f_max = 3;
sys.d_max = 2;
sys.to_life = 2; 

% x = [p_large, p_1, p_2, f_1, f_2, p_o, t_o, d_1, d_2]
grid_x = {1-0.5:n+0.5,...
          1-0.5:n+0.5,...
          1-0.5:n+0.5,...
          -2-0.5:sys.f_max+0.5,...
          -2-0.5:sys.f_max+0.5,...
          1-0.5:l+0.5,...
          0-0.5:sys.to_life+0.5,...
          0-0.5:sys.d_max+0.5,...
          0-0.5:sys.d_max+0.5}; 
sys.gpart_x = GridPartition(grid_x);

grid_u = {1-0.5:n+0.5,...
          1-0.5:n+0.5,...
          1-0.5:n+0.5}; 
sys.gpart_u = GridPartition(grid_u);


% load transitions
if strcmp(flag, 'load')
load('sys_tau_1to100.mat', 'sys_tau_1to100');
load('sys_tau_101to200.mat', 'sys_tau_101to200');
load('sys_tau_201toend.mat', 'sys_tau_201toend');
sys.tau = [sys_tau_1to100, sys_tau_101to200, sys_tau_201toend];
return; 
end

N = length(sys.gpart_x) + 1
M = length(sys.gpart_u);

% compute transitions

[I,J,K,S]=deal(cell(N-1,1));
    mycell = {};
for k = 1:N
    mycell{k} = 0;
end
progressStepSize = 1000;
ppm = ParforProgMon('Example', N, progressStepSize);

parfor (i = 1:N-1) % opts)
    if mod(i,progressStepSize)==0
        ppm.increment();
    end
    x_i = sys.gpart_x.idx2mid(i);
    p0 = x_i(1);
    p1 = x_i(2);
    p2 = x_i(3);
    f1 = x_i(4);
    f2 = x_i(5);
    po = x_i(6); 
    to = x_i(7);
    d1 = x_i(8);
    d2 = x_i(9);
    
    P0_plus = find(Map(p0,:));
    P1_plus = find(Map(p1,:));
    P2_plus = find(Map(p2,:));
    mytemp_p_k = sparse(N, M);
    for p0_plus = P0_plus
        for p1_plus = P1_plus
            for p2_plus = P2_plus
                k = sys.gpart_u.pnt2idx([p0_plus;p1_plus;p2_plus]);
                
                unsafe_ik = 0;
                if ~ismember(p0, sys.p_large_safe)    
                    % large car is not safe 
                    mytemp_p_k(N,k) = 1;
                else
                    if Map(p0,p0_plus) && Map(p1,p1_plus) && Map(p2,p2_plus) 
                        % valid move on the Map
                        
                        % dead time update
                        if ismember(p1,setdiff(1:n,sys.p_large_safe))
                            d1_plus = d1 - 1;
                        else
                            d1_plus = sys.d_max;
                        end
                        if ismember(p2,setdiff(1:n,sys.p_large_safe))
                            d2_plus = d2 - 1;
                        else
                            d2_plus = sys.d_max;
                        end
                        if d1_plus <= -0.5 || d2_plus <= -0.5 
                            mytemp_p_k(N,k) = 1;
                            unsafe_ik = 1;  
                        end

                        % fuel update
                        if ~unsafe_ik
                            if f1 >= 0.5 % enough fuel for nontrivial move
                                if p1 ~= p1_plus
                                    f1_plus = f1 - 1;
                                else
                                    f1_plus = f1;
                                end
                            else % zero fuel
                                if p1 ~= p1_plus
                                    mytemp_p_k(N,k) = 1;
                                    unsafe_ik = 1;
                                else
                                    f1_plus = f1 - 1;
                                end
                            end

                            if f2 >= 0.5
                                if p2 ~= p2_plus  
                                    f2_plus = f2 - 1;
                                else
                                    f2_plus = f2;
                                end
                            else
                                if p2 ~= p2_plus  
                                    mytemp_p_k(N,k) = 1;
                                    unsafe_ik = 1;
                                else
                                    f2_plus = f2 - 1;
                                end
                            end
                        end

                        if ~unsafe_ik
                            if f1_plus <= -2.5 || f2_plus <= -2.5
                                % over-draw fuel
                                    mytemp_p_k(N,k) = 1;
                                % unsafe_ik = 1;
                            else
                                % re-fuel from the large vehicle "v0"
                                if p1_plus == p0_plus
                                    f1_plus = sys.f_max;
                                end
                                if p2_plus == p0_plus
                                    f2_plus = sys.f_max;
                                end

                                if to >= -0.5 && (p1 == po || p2 == po)
                                    % capture now
                                    to_plus = sys.to_life; 
                                    po_plus = sys.p_target;
                                else
                                    % capture later
                                    to_plus = to - 1;
                                    po_plus = po;
                                end

                                if to_plus <= -0.5
                                    % target dies before being captured
                                    mytemp_p_k(N,k) = 1;
                                    % unsafe_ik = 1;
                                else
                                    my_i_plus = [];
                                    for p = po_plus
                                        i_plus = sys.gpart_x.pnt2idx([p0_plus; p1_plus; p2_plus; f1_plus; f2_plus; p; to_plus; d1_plus; d2_plus]);
                                        my_i_plus = [my_i_plus; i_plus];
                                        
                                    end
                                    ind = sub2ind(size(mytemp_p_k), i_plus, k*ones(size(i_plus)));
                                    mytemp_p_k(ind) = 1;
                                end
                            end
                        end
                    else
                        % invalid move on the Map
                        %h = h+1;
                        %mytemp(k,h) = N;
                        mytemp_p_k(N,k) = 1;
                        %tau{k}(i, N) = 1;
                    end
                end
            end
        end
    end
mycell{i} = mytemp_p_k;
[I{i},J{i},S{i}]=find(mytemp_p_k);
K{i}=repmat(i,size(S{i}));
end
ppm.delete();

I=cell2mat(I(:));
J=cell2mat(J(:));
K=cell2mat(K(:));
S=cell2mat(S(:));
tau=ndSparse.build([I,J,K],S);
tau(N,:,N) = 1;

sys.tau = tau;
